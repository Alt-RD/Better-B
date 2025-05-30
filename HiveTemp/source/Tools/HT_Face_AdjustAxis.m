%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University, AltRD-Emmanuel Ruffio
%  Author: emmanuel.ruffio@alt-rd.com
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
% Adjust B axis and change coordinates so that B alignement matches A alignement
% or custom vectors u,v provided
% It is useful to compute distance between nodes of different faces
% - A,B: two face structure
%   A can be dim 3x3 specifying a coordinate system. In that case, u will be
%     changed to match the first possible vector in that system. v will be set
%     so that u,v,norm are a direct coordinate system (norm is the vector pointing
%     outside the medium).
%   A can be dim 3x2 specifying two vectors u,v. In that case, coordinate system
%     of the face is changed so that u vectors match, v as well.
%   A can be dim 3x1 specifying only vector u. In that case, coordinate system
%     of the face is changed so that u vector match specified vector.
% - uvB: a two/four component row vector [0;1] usually used to select a subpart of face B
%   uvB will be transformed identically to face B. Multiple rows will be transformed
%   similarly. Set uvB to [] if not used
% - options:
% Returns
% 1) the new face B
% 2) the face selection vector modified the same way as face axis
% Changed:
% 15/06/2021: Add possiblity of specifying A dim3x1 only.
function [B uvB] = HT_Face_AdjustAxis(B, A, uvB, options)
  assert(nargin >= 2);

  if nargin < 3, uvB = []; endif;
  if nargin < 4, options = struct('verbose', true); endif

  assert(HT_CheckType(B, 'face'), 'Invalid face object B');
  assert(~isempty(B.axis), 'Invalid face object B. No axis is defined'); % COuld happen with an empty face
  assert(isempty(uvB) || (size(uvB, 2) == 2) || (size(uvB, 2) == 4));
  assert(isempty(uvB) || (all(uvB <= 1) && all(uvB>=0)));

  X=1; Y=2; Z=3;

  if ~isstruct(A)
    if all(size(A) == [3 3])
      lNewAxis = [];
      if abs(dot(B.norm, A(:,X)) - 1) < 1E-10 % Face B is normal to X and pointing towards X ?
        lNewAxis = A(:, [Y Z]);
      elseif abs(dot(B.norm, A(:,X)) + 1) < 1E-10 % Face B is normal to X and pointing towards -X ?
        lNewAxis = A(:, [Y -Z]);
      elseif abs(dot(B.norm, A(:,Y)) - 1) < 1E-10 % Face B is normal to Y and pointing towards Y ?
        lNewAxis = A(:, [X -Z]);
      elseif abs(dot(B.norm, A(:,Y)) + 1) < 1E-10 % Face B is normal to Y and pointing towards -Y ?
        lNewAxis = A(:, [X Z]);
      elseif abs(dot(B.norm, A(:,Z)) - 1) < 1E-10 % Face B is normal to Z and pointing towards Z ?
        lNewAxis = A(:, [X Y]);
      elseif abs(dot(B.norm, A(:,Z)) + 1) < 1E-10 % Face B is normal to Z and pointing towards -Z ?
        lNewAxis = A(:, [X -Y]);
      else
        error('Coordinate system or face coordinate system is invalid');
      endif

      A = HT_Face_Init(sprintf('custom_u=(%d,%d,%d)_v=(%d,%d,%d)', lNewAxis), ...
                         'axis', lNewAxis);
      clear lNewAxis;
    elseif all(size(A) == [3 2])
      A = HT_Face_Init(sprintf('custom_u=(%d,%d,%d)_v=(%d,%d,%d)', A), ...
                          'axis', A);
    elseif all(size(A) == [3 1])
      if     abs(abs(dot(B.axis(:,1), A)) - 1) < 1E-10 % Axis(U) of face B is aligned with A
        A = HT_Face_Init(sprintf('custom_u=(%.1f,%.1f,%.1f)', A), ...
                            'axis', [A B.axis(:,2)]);
      elseif abs(abs(dot(B.axis(:,2), A)) - 1) < 1E-10 % Axis(V) of face B is aligned with A
        A = HT_Face_Init(sprintf('custom_u=(%.1f,%.1f,%.1f)', A), ...
                            'axis', [A B.axis(:,1)]);
      else
        error(sprintf('Vector u=(%.1f,%.1f,%.1f) specified but it does not match face axis', A));
      endif
    else
      error('Invalid parameter <A> specifying vectors u,v');
    endif
  endif

  assert(~isempty(A.axis), 'Invalid face object A. No axis is defined'); % COuld happen with an empty face

  % Compute the normal vector of the face
  zB = cross(B.axis(:,1), B.axis(:,2));
  zA = cross(A.axis(:,1), A.axis(:,2));

  zDot = dot(zA, zB);

  if abs(zDot) < 1E-14
    % If faces are not in the same plane, the user must specified the axis explicitly
    error(sprintf('Faces <%s.%s> and <%s.%s> are not in the same plane', A.model, A.name, B.model, B.name));
  else %if they are in the same plane
    % If uA == uB are not aligned, uB, and vB are flip
    if abs(dot(A.axis(:,1), B.axis(:,1))) < 1E-10
      B.size = flip(B.size);
      B.axis = fliplr(B.axis);
      B.n = flip(B.n);
      B.pos = fliplr(B.pos);

      if isempty(B.vertex) % Rectangles based face ?
        if ~isempty(B.dims), B.dims = B.dims(:,[3 4 1 2]); endif %[B.dims(:,[3 4]) B.dims(:,[1 2])]; endif
      else % General based face ?
        B.vertex = B.vertex(:,[2 1]);
      endif

      if size(uvB,2) == 2, uvB = fliplr(uvB); endif;
      if size(uvB,2) == 4, uvB = uvB([3 4 1 2]); endif
    endif

    % If uA and uB are reversed, uB is reversed
    if dot(A.axis(:,1), B.axis(:,1)) < 0
      B.globalPosition += B.size(1) * B.axis(:,1);
      B.axis(:,1) = -B.axis(:,1);

      if ~isempty(B.pos), B.pos(:,1) = 1 - B.pos(:,1); endif
      if isempty(B.vertex) % Rectangles based face ?
        if ~isempty(B.dims), B.dims(:,[1 2]) = 1 - B.dims(:,[2 1]); endif
      else
        B.vertex(:,1) = 1 - B.vertex(:,1);
      endif

      if size(uvB,2) == 2, uvB(:,1) = 1 - uvB(:,1); endif;
      if size(uvB,2) == 4, uvB(:,[1 2]) = 1 - uvB(:,[2 1]); endif
    endif

    % If vA and vB are reversed, vB is reversed
    if dot(A.axis(:,2), B.axis(:,2)) < 0
      B.globalPosition += B.size(2) * B.axis(:,2);
      B.axis(:,2) = -B.axis(:,2);

      if ~isempty(B.pos), B.pos(:,2) = 1 - B.pos(:,2); endif
      if isempty(B.vertex) % Rectangles based face ?
        if ~isempty(B.dims), B.dims(:,[3 4]) = 1 - B.dims(:,[4 3]); endif
      else
        B.vertex(:,2) = 1 - B.vertex(:,2);
      endif

      if size(uvB,2) == 2, uvB(:,2) = 1 - uvB(:,2); endif;
      if size(uvB,2) == 4, uvB(:,[3 4]) = 1 - uvB(:,[4 3]); endif
    endif
  endif


endfunction
