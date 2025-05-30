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
% Returns the flux factor, i.e. the fraction of incoming heat radiation
% absorbed by hive faces.
%
% Input arguments
% cosBeta = (dim Nt x Nf) cosine of angle between face normal and sun position
%           Nt refers to time sample and Nf to face index
% epsf = (dim 1 or dim Nf) fonction(angle) or real value: emissivity with respect to incoming
%           radiation angle
%
% Retour:
% A (dim Nt x Nf): Fraction of heat flow absorbed by each face and at every time step
%                  Dimension of A is similar to input argument cosBeta.
function [A] = HT_GetSolarFluxFactor(cosBeta, epsf)

  nf = size(cosBeta, 2);

  % Checks that emissivity argument is valid
  assert(is_function_handle(epsf) || isa(epsf, 'numeric') || iscell(epsf), 'Invalid argument <epsf>');
  % Size of vector epsf must be 1 or size(cosBeta, 2) (number of faces)
  if isa(epsf, 'numeric')
    assert(any(numel(epsf) == [1, nf]), 'Invalid vector size <epsf>');
    if (numel(epsf) == 1), epsf = repmat(epsf, nf, 1); endif;

    lValue = epsf;
    epsf = cell(nf, 1);
    for i=1:nf
      epsf{i} = @(x) lValue(i);
    endfor
  elseif is_function_handle(epsf)
    lValue = epsf;
    epsf = cell(nf, 1);
    for i=1:nf
      epsf{i} = lValue;
    endfor
  else
    assert(numel(epsf) == nf, 'Invalid cells size <epsf>');
    assert(all(cellfun(@is_function_handle, epsf)), 'Invalid cells content <epsf>');
  endif

  lPositiveAngle = (cosBeta > 0);
  A = cosBeta;
  A(~lPositiveAngle) = 0;

  for i=1:nf
    A(lPositiveAngle(:,i), i) *= epsf{i}( acos(cosBeta(lPositiveAngle(:,i), i)) );
  endfor

end


%  lCosAngleZ = cos(solGamma)*cos(lSolBeta);
%  lCosAngleXp = sin(solGamma)*cos(alpha) - cos(solGamma)*sin(lSolBeta)*sin(alpha);
%  lCosAngleXm = sin(solGamma)*cos(pi+alpha) - cos(solGamma)*sin(lSolBeta)*sin(pi+alpha);
%  lCosAngleYp = sin(solGamma)*cos(pi/2+alpha) - cos(solGamma)*sin(lSolBeta)*sin(pi/2+alpha);
%  lCosAngleYm = sin(solGamma)*cos(alpha-pi/2) - cos(solGamma)*sin(lSolBeta)*sin(alpha-pi/2);
%
%  C = [lCosAngleZ; lCosAngleXp; lCosAngleXm; lCosAngleYp; lCosAngleYm];
%  A = acos(C);
