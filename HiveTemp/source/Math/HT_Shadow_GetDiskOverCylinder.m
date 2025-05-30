%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University, AltRD-Emmanuel Ruffio
%
%  HiveTemp is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  HiveTemp is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with HiveTemp.  If not, see <https://www.gnu.org/licenses/>
% ========================================================================
% Returns the fraction of a cylinder that is hidden below a disk
% This function is used to compute the shadow of a circular stone over
% a trunc hive
%
% Input arguments
% diskRadius = (dim 1 x 1) radius of the upper disk
% cRadius = (dim 1 x 1) radius of the cylinder below the disk
% heightVec = (dim N x 1) [h0, h1, ... hN] start and end position of each cylinder
%                       part. The disk is assumed to be at z=0, and the cylinder
%                       looking backward (to Z-), so h1 > h0
% cTheta = (dim L x 1) [theta0, theta1, ... thetaN] subdivision of the cylinder
%                       along the theta direction (default is [0, 2pi])
% posMat = (dim 2 x M)  direction of the light (radian).
%                       [azimuth; elevation] angles.
%                       Azimuth 0 means angle are not changed. It means azimuth is referenced
%                       to the first axis used to define the cylinder and disk geometry
%
% Options:
% .sort = ['theta'] or 'height': defines how data are packed into matrices.
%                      'theta'  => [(theta1,height1) ; (theta2,height1); (thetaN;height1) ...]
%                      'height' => [(theta1,height1) ; (theta1,height2); (theta1;heightN) ...]
%
% Retour:
% A (dim (NxL) x M): Illuminated part of each cylinder part, ordered by L first, then by N
% V (dim (NxL) x M): Solar heat flux absorbed by each cylinder part, ordered by L (theta) first, then by N (height)
% S (dim (NxL) x 1): Node area
function [A V S] = HT_Shadow_GetDiskOverCylinder(diskRadius, cRadius, heightVec, thetaVec, posMat, varargin)
  HT_ImportConstants();

  lTicId = tic();

  % Check input parameters
  assert(nargin >= 5,                                 'Shadow_GetDiskOverCylinder::Missing input arguments');

  assert(isscalar(diskRadius) && (diskRadius >= 0),   'Shadow_GetDiskOverCylinder::Invalid disk radius');
  assert(isscalar(cRadius) && (cRadius >= 0),         'Shadow_GetDiskOverCylinder::Invalid cylinder radius');
  assert(diskRadius >= cRadius,                       'Shadow_GetDiskOverCylinder::Invalid roof radius');
  assert(~isempty(heightVec) && all(heightVec >= 0),  'Shadow_GetDiskOverCylinder::Invalid cylinder height');
  % If node are too large (> pi/2), some issues could appear in the algorithm
  assert(isempty(thetaVec) || (isnumeric(thetaVec) && all(diff(thetaVec) < 2*pi/5)), 'Shadow_GetDiskOverCylinder::Invalid cylinder theta subdivision');
  assert(all(posMat(2,:) <= pi/2),                    'Shadow_GetDiskOverCylinder::Invalid elevation angle. Must be lower than pi/2');
  assert(rows(posMat) == 2,                           'Shadow_GetDiskOverCylinder::Invalid sun position matrix');

  if isempty(thetaVec), thetaVec = [-pi pi]; endif;
  if ~iscolumn(thetaVec), thetaVec = thetaVec'; endif;
  if ~iscolumn(heightVec), heightVec = heightVec'; endif;

  % If <heightVec> is scalar, implicitly preprend a 0
  if isscalar(heightVec), heightVec = [0; heightVec]; endif;
  heightVec -= heightVec(end);

  lParams = struct();
  lOptions = struct();
  lParamEndIndex = numel(varargin);

  if numel(varargin) > 1
    if isstruct(varargin{end}) && isstruct(varargin{end-1})
      lParams = varargin{end-1};
      lOptions = varargin{end};
      lParamEndIndex -= 2;
    elseif isstruct(varargin{end})
      lParams = varargin{end};
      lParamEndIndex -= 1;
    endif
  endif

  lOptions = HT_CheckField(lOptions, 'verbose',        true,     @(v) islogical(v));
  lOptions = HT_CheckField(lOptions, 'verboselevel',   true,     { @(v) isscalar(v) && (v >= 0) && (round(v) == v) } );
  lOptions = HT_CheckField(lOptions, 'skipFieldCheck', false,    @(v) islogical(v));
  lOptions = HT_CheckField(lOptions, 'sort',           'theta',  @(v) any(strcmpi(v, {'theta', 'height'})));
  lOptions = HT_CheckField(lOptions, 'useCache',       true,      @(v) islogical(v));

  if lParamEndIndex > 0
    prop = varargin(1:2:lParamEndIndex);
    values = varargin(2:2:lParamEndIndex);

    assert(numel(prop) == numel(values), 'Invalid input parameters');

    % Add to the <lParams> structure all parameters specified
    for i=1:numel(prop)
      if any(strcmpi(prop{i}, {'sort'}))
        lOptions = setfield(lOptions, 'sort', values{i});
      else
      lParams = setfield(lParams, prop{i}, values{i});
      endif
    endfor
  endif

  % Check parameter structure fields. It is desirable to check that no silly errors
  % are made on the parameter name
  if ~lOptions.skipFieldCheck
    lValidParamList = { 'power', 'epsType', 'deltaDisk' };
    lParamList = fieldnames(lParams);
    lValidParam = cellfun(@(v) any(strcmp(v, lValidParamList)), lParamList);
    assert(all(lValidParam), sprintf('Invalid parameter name detected <%s>', strjoin(lParamList(~lValidParam), ',')));
    clear lParamList lValidParam lValidParamList;
  endif

  lParams = HT_CheckField(lParams, 'power',             0,                 {@(v) isscalar(v) && (v >= 0) });
  lParams = HT_CheckField(lParams, 'epsType',           'constant',        {@(v) any(strcmpi(v, {'constant', 'cosine', 'power'})) });
  lParams = HT_CheckField(lParams, 'deltaDisk',         [0, 0],            {@(v) numel(v)==2 });

  assert(norm(lParams.deltaDisk, 2) + cRadius <= diskRadius, 'The shift of the roof is too large');

  lReload = false;

  % Attempt to reload data saved in cache file
  lCacheFileName = 'HT_Shadow_GetDiskOverCylinder_Cache.mat';
  if lOptions.useCache && exist(lCacheFileName, 'file')
    D = load(lCacheFileName);
    lReload = isfield(D, 'lOptions') && isfield(D, 'lParams');
    lReload = lReload && isequal(D.lOptions, lOptions) && isequal(D.lParams, lParams);
    lReload = lReload && isfield(D, 'A') && isfield(D, 'S') && isfield(D, 'V');

    if lReload
      if lOptions.verbose
        disp('HT_Shadow_GetDiskOverCylinder::Data have not changed. Reload from cache file.');
      endif

      A = D.A;
      V = D.V;
      S = D.S;
    endif
  endif

  if ~lReload
    lSortTheta = strcmpi(lOptions.sort, 'theta');

  nHeight = numel(heightVec);
  nTheta = numel(thetaVec);
  nAngle = columns(posMat);

  heightDim = [heightVec(1:(end-1)), heightVec(2:end)] / cRadius; % Normalized height

  A = zeros((nHeight-1)*(nTheta-1), nAngle);
  V = zeros((nHeight-1)*(nTheta-1), nAngle);

  dDisk = lParams.deltaDisk/cRadius; % Decalage normalise du toit

  cSpace = heightVec(2:end) - heightVec(1:(end-1));
  alpha = pi/2/(pi/2)^lParams.power;

    S = cRadius * (thetaVec(2:end)-thetaVec(1:(end-1))) * (heightVec(2:end) - heightVec(1:(end-1)))';

  % Initialize the emissivity(theta) function, that handles the dependency
  % of the emissivity to the light angle
  if strcmpi(lParams.epsType, 'constant')
    epsFunc = @(beta) 1;
  elseif strcmpi(lParams.epsType, 'power')
    epsFunc = @(beta) (1-(lParams.power > 0)*(beta/(pi/2)).^lParams.power);
  else
    epsFunc = @(beta) cos(beta.^lParams.power * (lParams.power > 0) * alpha);
  endif

  for i=1:nAngle
    % First, rotate the cylinder mesh so that the axis match the sun position
    % The only relevant angle is then the elevation angle
    lAzimuthAngle = posMat(1,i);
    lElevationAngle = posMat(2,i);
    lSunDx = cos(lElevationAngle);

    % If the sun is below the horizon, the cylinder is assumed to be shadowed
    if lElevationAngle < 0
      continue;
    endif

    lCurrentThetaVec = thetaVec - lAzimuthAngle;
    % Fix angle reference so that all angles lie in [-pi;pi];
    lCurrentThetaVec = mod(lCurrentThetaVec, 2*pi);
      lFaceThetaDim = [lCurrentThetaVec(1:(end-1)), lCurrentThetaVec(2:end)];
      lFaceThetaDim(lFaceThetaDim(:,1) >= lFaceThetaDim(:,2), 2) += 2*pi;
      lFaceThetaDim((lFaceThetaDim(:,2) > 3/2*pi),:) -= 2*pi;
      lFaceThetaDim = clip(lFaceThetaDim, [-pi/2, pi/2]);

    % Only the faces oriented towards the sun will be illuminated
    % thetaVec in [-pi/2, pi/2]
      lFaceIlluminatedFlag = (lFaceThetaDim(:,2) != lFaceThetaDim(:,1));

    % Theta is the angle on the cylinder, over which the integration is done
    % The normal vector is aligned with this direction ( cos(theta), sin(theta), 0 )
    % since the hive is horizontal
    % Voir rapport ruche Tronc
    % J = integrate  dot(phi,n) * epsilon(beta) * r* dtheta
    % J = integrate  dx.cos(theta) * epsilon(beta) * r* dtheta
    % Avec dx = cos(lElevationAngle)
    % Constant are not inserted in the function and must be integrated after the ingration
    % => (z(2)-z(1));
    irrFunc = @(theta) cRadius*lSunDx*cos(theta) .* epsFunc(acos(lSunDx*cos(theta)));

    % For each angle, the position of the shadow line is computed
    fHeight = @(theta) tan(lElevationAngle) * (cos(theta) - sqrt((diskRadius/cRadius)^2 - (sin(theta) - lParams.deltaDisk(2)/cRadius).^2) - lParams.deltaDisk(1)/cRadius);
    lShadowHeight = fHeight(lFaceThetaDim);

    % Determine l'angle d'inversion de la courbe de hauteur (dérivée dz/dtheta = 0)
    df=@(t) (cos(t).*(sin(t) - dDisk(2)) ./ sqrt((diskRadius/cRadius)^2 - (sin(t) - dDisk(2))^2) - sin(t)).^2;
    [lThetaZero, fval, info] = fminbnd(df, -pi/2, pi/2, ...
                                           optimset('TolX', 1E-6, ...
                                                    'TolFun', 1E-10));
    assert(info == 1, 'The algorithm did not converge');
    assert(fval < 1E-10, 'The solution found is not a zero. The position shift of the roof may be too large.');
    lThetaZero = round(lThetaZero/(2*pi*HT_VAR_EPSILON_U))*2*pi*HT_VAR_EPSILON_U;
    lThetaZeroHeight = fHeight(lThetaZero);

    for k=1:(nTheta-1) % nTheta is the number of vertices, not the number of face
      if ~lFaceIlluminatedFlag(k), continue; endif

      % If the maximum of the shadow line occurs in this set [theta0, theta1]
      % the maximum height of the shadow line is not anymore at the edge but
      % is somewhere in the middle of the theta interval.
      if (lFaceThetaDim(k,1) < lThetaZero) && (lThetaZero < lFaceThetaDim(k,2))
        lShadowHeightMin = min(lShadowHeight(k,:));
          lShadowHeightMax = max([lShadowHeight(k,:) , lThetaZeroHeight]);
      else
        lShadowHeightMin = min(lShadowHeight(k,:));
        lShadowHeightMax = max(lShadowHeight(k,:));
      endif

      lFullyIlluminated = (heightDim(:,2) <= lShadowHeightMin);
      lFullyOccluded =    (heightDim(:,1) >= lShadowHeightMax);

      [lFullIrradiance, ier] = quad(irrFunc, lFaceThetaDim(k,1), lFaceThetaDim(k,2));
      assert(ier == 0, sprintf('Could not integrate irradiance over the interval [%.2f;%.2f]', lFaceThetaDim(k,1), lFaceThetaDim(k,2)));

      lIrradiance = NA(nHeight-1,1);
      lIrradiance( lFullyIlluminated ) = lFullIrradiance * (heightDim(lFullyIlluminated,2) - heightDim(lFullyIlluminated,1));
      lIrradiance( lFullyOccluded ) = 0;

      % Ratio of shadow part area over node area for each nodes
      lIlluminatedRatio = NA(nHeight-1,1);
      lIlluminatedRatio( lFullyOccluded ) = 0;
        lIlluminatedRatio( lFullyIlluminated ) = cRadius^2*((lFaceThetaDim(k,2) - lFaceThetaDim(k,1))*(heightDim(lFullyIlluminated,2) - heightDim(lFullyIlluminated,1))) ./ ...
                                                        S(k,lFullyIlluminated)';

      for m=find(isna(lIlluminatedRatio))'
        [lIlluminatedRatio(m) lIrradiance(m)] = Int_IntegrateShadow(tan(lElevationAngle), ...
                                                                    diskRadius/cRadius, ...
                                                                    dDisk, ...
                                                                    lFaceThetaDim(k,:)', ...
                                                                    heightDim(m,:), ...
                                                                    lThetaZero, ...
                                                                    irrFunc);
      endfor

        if lSortTheta
      t = k + (0:(nHeight-2))*(nTheta-1);
        else
          t = (k-1)*(nHeight-1) + (1:(nHeight-1));
        endif

      A(t, i) = lIlluminatedRatio;
      V(t, i) = lIrradiance * cRadius;    % cRadius converts normalized height to real height

    endfor
  endfor

    % Surface matrix S is ordered by theta first (by default)
    if ~lSortTheta
      S = S';
    endif

    S = S(:); % Convert surface matrix (theta x height) into column vector

    if lOptions.useCache
      save('-binary', lCacheFileName, 'lOptions', 'lParams', 'A', 'S', 'V');
    endif
  endif % If ~lReload

  if lOptions.verbose, disp(sprintf('Computing shadow on cylinder done in %.2f ms', toc(lTicId))); endif; % Start timer (tic)

end

% Input parameters:
% hVec: [z0 z1] height boundaries for integration
% thetaVec: [theta0 theta1] angle boundaries for integration
% zVec : [zmin zmax] low and highest height of the projected shadow
% zFunction: @(theta) function which returns the height of the projected shadow for a given angle
% epsilonThetaFunc: @(theta) function which returns the emissivity for given angle
function V = Int_IntegrateSun(hVec, thetaVec, zVec, zFunction, epsilonThetaFunc)
  z1 = zFunction(thetaVec(2));
  V = quadv(@(theta) cos(theta) .* epsilonThetaFunc(theta) .* (z1 - zFunction(theta)) , thetaVec(1), thetaVec(2), ...
            1E-6, false);

  if thetaVec(1) > 0 % Some part of the node is fully illuminated
    V += (min(hVec(2), zVec(2))-hVec(1)) * quadv(@(x) cos(x).*epsilonThetaFunc( x ), 0, thetaVec(1));
  endif
endfunction

% Objective function used to compute the angle theta that correponds to a height z.
% Input parameters:
% s: s=-tan(elevation angle)
% a: ratio of disk radius R and cylinder radius r => a=R/r
% h: normalized height => h=z/r;
function [J M] = Int_ComputeShadowAngleSolution(s, a, h, theta)
  sinTheta = sin(theta);
  sinThetaSqr = sinTheta.^2;
  cosTheta = cos(theta);

  J = (s./h) .* (cosTheta - sqrt(a^2-sinThetaSqr)) - 1;

  if nargout > 1
    M = diag((s./h).*sinTheta.*(cosTheta./sqrt(a^2-sinThetaSqr)-1));
  endif
endfunction

% Integrate the shadow over theta=[thetaSet(1), thetaSet(2)] and z=[zSet(1), zSet(2)];
% If the theta boundaries lead to z values outside zSet, theta
% is adjusted
% Input parameters:
% -> s = ez/ex
% -> r = disk radius over cylinder radius
% -> d = [dim 2x1] delta position normalized by cylinder radius
% -> thetaSet = [dim 2x1] (radian) the angle set of integration
% -> zSet = [dim 2x1] (-) normalized height by cylinder radius. Ascending order (zSet(2) > zSet(1))
% -> thetaZero = [double] in [-pi/2;pi/2] theta that nullifies the derivative dz/dtheta
% -> irrFunc = @(theta) a function that gives the solar irradiance with respect to angle theta
%                       The surface emissivity is integrated in this function.
% Output:
% -> V = the ratio of shadowed area of the cylinder part
function [V J] = Int_IntegrateShadow(s, r, d, thetaSet, zSet, thetaZero, irrFunc)
  assert(iscolumn(thetaSet));

  f=@(theta) s*(cos(theta) - sqrt(r^2 - (sin(theta) - d(2)).^2) - d(1)) - zSet(1);

  % irrFunc does not take into account the height of integration z(theta)
  irrFuncZ=@(theta) irrFunc(theta) .* f(theta);

  dz = zSet(2)-zSet(1);
  S = (thetaSet(2)-thetaSet(1))*dz;
  V = 0; % Will count the area of the shadowed part
  J = 0; % Will count the energy absorbed by the area.

  % if thetaZero lies in the interval thetaSet, thetaSet must be splitted in
  % two parts.
  % Voir page 20 cahier de brouillon
  if (thetaZero > thetaSet(1)) && (thetaZero < thetaSet(2))
    thetaSet = [thetaSet(1) thetaZero;
                thetaZero, thetaSet(2)];
  endif

  nIntegral = columns(thetaSet);
  for i=1:nIntegral
    % Now, we adapt theta so that z(theta) lies in the specified interval zSet
    originThetaSet = thetaSet(:,i);
    clampThetaSet = thetaSet(:,i);
    zBoundary = f(clampThetaSet);

    if all(zBoundary <= 0) % Shadow is lower than the z set.
      V += (originThetaSet(2)-originThetaSet(1))*dz;    % The area is entirely covered by the shadow
    elseif all(zBoundary >=dz) % Node is fully illuminated
      V += 0;
      [Jpart, ier] = quad(irrFunc, originThetaSet(1), originThetaSet(2));
      assert(ier == 0, 'Could not integrate solar irradiance over the cylinder');
      J += Jpart*dz;
    else
      if zBoundary(1) < 0       % If the shadow line is lower than the left corner
        clampThetaSet(1) = Int_GetThetaFromZ(s,r,d,zSet(1), clampThetaSet);
        V = (clampThetaSet(1) - originThetaSet(1))*dz;

      elseif zBoundary(2) < 0   % If the shadow line is lower than the right corner
        clampThetaSet(2) = Int_GetThetaFromZ(s,r,d,zSet(1), clampThetaSet);
        V = (originThetaSet(2) - clampThetaSet(2))*dz;

      elseif zBoundary(1) > dz     % If the shadow line is higher than the left corner
        clampThetaSet(1) = Int_GetThetaFromZ(s,r,d,zSet(2), clampThetaSet);

        [Jpart, ier] = quad(irrFunc, originThetaSet(1), clampThetaSet(1));
        assert(ier == 0, 'Could not integrate solar irradiance over the cylinder');
        J += Jpart*dz;
      elseif zBoundary(2) > dz     % If the shadow line is higher than the right corner
        clampThetaSet(2) = Int_GetThetaFromZ(s,r,d,zSet(2), clampThetaSet);

        [Jpart, ier] = quad(irrFunc, clampThetaSet(2), originThetaSet(2));
        assert(ier == 0, 'Could not integrate solar irradiance over the cylinder');
        J += Jpart*dz;
      endif

      [I err] = quad(f, clampThetaSet(1), clampThetaSet(2));
      I -= dz*(clampThetaSet(2) - clampThetaSet(1));
      assert(err == 0, 'Could not integrate solar irradiance over the cylinder');
  ##    I -= zSet(1)*(clampThetaSet(2)-clampThetaSet(1));
      V -= I;

      [I err] = quad(irrFuncZ, clampThetaSet(1), clampThetaSet(2));
      assert(err == 0, 'Could not integrate solar irradiance over the cylinder');
      J += I;
    endif
  endfor

  V = 1 - V/S;

endfunction


% Input parameters:
% -> s = ez/ex
% -> r = disk radius over cylinder radius
% -> d = [dim 2x1] delta position normalized by cylinder radius
function thetaSol = Int_GetThetaFromZ(s, r, d, zValue, thetaSet)
##  fz = @(t) (s*(cos(t) - sqrt(r^2-(sin(t)-d(2)).^2) - d(1)) - zValue).^2;
##
##  [thetaSol, fval, info, output] = fminbnd(fz, thetaSet(1), thetaSet(2), ...
##                            optimset( 'Jacobian', 'off', ...
##                                      'MaxFunEvals', 100, ...
##                                      'TolX', 1E-6, ...
##                                      'TolFun', 1E-10, ...
##                                      'MaxIter', 20, 'display', 'off'))
##  fz = @(t) (s*(cos(t) - sqrt(r^2-(sin(t)-d(2)).^2) - d(1)) - zValue).^2;

% I did not see fzero support from search boundaries
  fz = @(t) s*(cos(t) - sqrt(r^2-(sin(t)-d(2)).^2) - d(1)) - zValue;
  [thetaSol, fval, info, output] = fzero(fz, thetaSet, ...
                            optimset( 'MaxFunEvals', 100, ...
                                      'TolX', 1E-6, ...
                                      'TolFun', 1E-10, ...
                                      'MaxIter', 20, 'display', 'off'));

   assert(info == 1, 'Failed to converge properly');
endfunction
