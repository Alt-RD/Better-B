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
% Return the cosine of angles between hive faces and the sun.
%
% Input arguments
% 1) phi   = zenith angle of the sun (dim N x 1) [in radian]
% 2) theta = azimuth angle of the sun (relative to geo. North and clockwise)
%         (dim N x 1) [in radian]
% 3) vec   = specifies the orientation of the hive. options specifies if it
%         is a list of normal vectors or a 2-dim vector (azimuth, inclination)
%         Normal vectors are always pointing outside the hive.
%
%         - if options is '':  matrix (Nn x 3) containing normal vectors of the 
%                            beehive faces
%         - if options.angle is true: matrix (1 x 2) = [alpha, gamma]
%               with alpha = azimuth of the hive (0 is north, pi/2 is east)
%                    gamma = inclination of the hive (zenith angle - pi/2)
%                            0 means the hive is horizontal
%               both angles refers to the hive normal vector of entrance face.
% 4) options   = struct('paramname', paramvalue,...) containing options
%   - 'angle' -> = true/false, specify if vec is a 2-elt vectors (azimuth, inclination)
%                              or normal vectors (angle is false)
%
% Output arguments:
% 1) cosBeta (dim N x Nn): cosine of angles between solar positions(i in [1:N]) and face normals(k in [1:Nn])
%      If options.angle is true: Nn == 6 since 6 normals are generated
%                                in the following order
%                                 1) X- hive (opposite to entrance)  (ex: north)
%                                 2) X+ hive entrance                (ex: south)
%                                 3) Y- hive side                    (ex: east)
%                                 4) Y+ hive side (clockwise)        (ex: west)
%                                 5) Z- hive bottom side             (ex: towards ground)
%                                 6) Z+ hive top side                (ex: towards zenith)
%      else    Nn is unchanged and equal to size(vec,1)
% 2) vec: matrix (dim Nn, 3) containing normal vectors
function [cosBeta vec] = HT_GetHiveSolarAngles(phi, theta, vec, options)
  
  assert(nargin >= 3, 'Missing parameters');
  assert(numel(phi) == numel(theta), 'Arg<phi> and Arg<theta> must have the same size');
  
  if nargin < 4, options = struct(); endif
  if ~isfield(options, 'angle'), options.angle = false; endif;
  
  phi = phi(:);
  theta = theta(:);
  N = numel(phi);
  
  if options.angle
    % In that case, vec contains two angles (azimuth and inclination) of the
	  % hive. 
    assert(all(size(vec) == [1 2]));
    
    alpha = vec(1);
    gamma = vec(2);
    
    % 6 faces and 3 components of each vector
    angle = [ alpha+pi, -gamma; ... X-  
              alpha, gamma; ...     X+
              alpha-pi/2, 0; ...    Y-
              alpha+pi/2, 0; ...    Y+
              alpha+pi, pi/2-gamma; ... Z- 
              alpha, gamma-pi/2]; ...Z+
              
    Nn = 6;
    cosBeta = (sin(phi) * cos(angle(:,2)')) .* cos(repmat(theta, 1, Nn) - repmat(angle(:,1)', N, 1)) - cos(phi) * sin(angle(:,2)');
    %beta = acos(cosBeta);

    if nargout > 1
      vec = [cos(angle(:,2)) .* cos(angle(:,1)),     cos(angle(:,2)) .* sin(angle(:,1)),    -sin(angle(:,2))];
    endif
  else % Normal vectors ?
    assert(size(vec,2) == 3);
    
    solarDir = [sin(phi) .* cos(theta),    sin(phi).*sin(theta),    cos(phi)];
    cosBeta = solarDir * vec';    % Dot product
    %beta = acos(cosBeta);
  endif;
end

