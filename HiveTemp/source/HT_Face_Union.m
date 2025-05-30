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

function F = HT_Face_Union(name, Fs, varargin)
  assert(nargin >= 2, 'There must be at least 1 input arguments');

  prop = varargin(1:2:end);
  values = varargin(2:2:end);

  assert(numel(prop) == numel(values), 'Invalid input parameters');

  % Default parameter values
  lMergeMesh = false;
  lAllowIntersect = false;   % Faces could overlap
  lMustBeCoplanar = false;

  for i=1:numel(prop)
    if strcmpi(prop{i}, "mesh")
      assert(islogical(values{i}), 'Invalid parameter <mergeMesh>');
      lMergeMesh = values{i};
    elseif strcmpi(prop{i}, "coplanar")
      assert(islogical(values{i}), 'Invalid parameter <coplanar>');
      lMustBeCoplanar = values{i};
    elseif strcmpi(prop{i}, "allowIntersect")
      assert(islogical(values{i}), 'Invalid parameter <allowIntersect>');
      lAllowIntersect = values{i};
    else
      error(sprintf('Invalid input parameter <%s>', prop{i}));
    endif
  endfor

  if isempty(Fs)
    F = {};
    return;
  endif

  if ~isstruct(Fs) && ~iscell(Fs)
    error('Invalid array of faces');
  endif

  % Check input faces
  assert(all(arrayfun(@(v) HT_CheckType(v, 'face'), Fs)), 'Invalid input faces');

  % Check orientation
  F = Int_GetObject(Fs(1)); % GetObject returns the object stored in a cell (1x1)
  assert(all(arrayfun(@(v) norm(cross(F.norm, Int_GetObject(v).norm), Inf) < 1E-10, Fs)), 'Faces are not oriented in the same direction');

  % Check coplanarity
  if lMustBeCoplanar
    global HT_VAR_EPSILON_POS

    lPlanePos = dot(F.globalPosition, F.norm);
    assert(all(arrayfun(@(v) abs(dot(Int_GetObject(v).globalPosition, F.norm) - lPlanePos) < HT_VAR_EPSILON_POS, Fs(2:end))), 'Faces are not coplanar');
    clear lPlanePos;
  endif

  % Build a Rect2 object for each 3d face
  if iscell(Fs)
    lRectList = cellfun(@(v) HT_Rect2.from3dSpace(  v.globalPosition, ...
                                                    v.axis,
                                                    v.size,
                                                    [F.axis F.norm]), ...
                          Fs, 'UniformOutput', false);
  else
    lRectList = arrayfun(@(v) HT_Rect2.from3dSpace( v.globalPosition, ...
                                                    v.axis,
                                                    v.size,
                                                    [F.axis F.norm]), ...
                          Fs, 'UniformOutput', false);
  endif

  % Now check that the union of faces can be done
  lUnion = HT_Rect2.fromUnion(lRectList, 'allowEmptySpace', false, 'allowIntersect', lAllowIntersect);
  assert(~lUnion.isEmpty(), 'Face_Union::Could not build the union of faces');

  if ~lMergeMesh % only merge the geometry
    F = HT_Face_Init(name,  'size', lUnion.size(), ...
                            'axis', F.axis, ...
                            'norm', F.norm, ...
                            'globalPosition', lUnion.to3dSpace(lUnion.origin(), [F.axis F.norm], F.globalPosition));
  else
    for i=2:numel(Fs)
      F = Int_GetUnion(F, Int_GetObject(Fs(i)), lAllowIntersect);
    endfor
  endif
endfunction

% Return an object. If obj is a structure, it does nothing.
% If a cell, it returns the content.
function obj = Int_GetObject(obj)
  if iscell(obj)
    obj = obj{1};
  endif
endfunction

function F1 = Int_GetUnion(F1, F2, lAllowIntersect)
  global HT_VAR_EPSILON_POS
  global HT_VAR_EPSILON_U

  lUMerge = abs(dot(F2.globalPosition - F1.globalPosition, F2.axis(:,2))) < HT_VAR_EPSILON_POS;
  lVMerge = abs(dot(F2.globalPosition - F1.globalPosition, F1.axis(:,1))) < HT_VAR_EPSILON_POS;

  % Find the coordinate of F2 vertices in F1 axis
  F2 = HT_Face_AdjustAxis(F2, F1.axis);
  lDims1 = F1.dims .* F1.size([1 1 2 2])';
  lDims2 = F2.dims .* F2.size([1 1 2 2])';

  lF2Offset = F2.axis' * (F2.globalPosition - F1.globalPosition);
  lDims2 += lF2Offset([1 1 2 2])';
  lDims = round([lDims1; lDims2] ./ HT_VAR_EPSILON_POS) * HT_VAR_EPSILON_POS;

  ##  lDims2(:, 1:2) += dot(F2.globalPosition - F1.globalPosition, F2.axis(:,1));
##  lDims2(:, 3:4) += dot(F2.globalPosition - F1.globalPosition, F2.axis(:,2));

  lPos1 = F1.pos .* F1.size';
  lPos2 = F2.pos .* F2.size';
  lPos2 += lF2Offset';
  lPos = round([lPos1; lPos2] ./ HT_VAR_EPSILON_POS) * HT_VAR_EPSILON_POS;

  % If some dims are negative, it means the globalPosition must be moved
  lDimsMin = min(lDims(:,[1 3]));
  lDims -= lDimsMin([1 1 2 2]);     % Change the origin of Dims so that all values start at 0
  lGlobalPosition = F1.globalPosition + F1.axis * lDimsMin';
  lSize = max(lDims(:,[2 4]), [], 1)';    % Compute new face size
  lDims ./= lSize([1 1 2 2])';         % Normalized node positions

  if lUMerge
    assert(all(abs(HT_Face_GetAbsoluteData(F1, 'vgrid') - HT_Face_GetAbsoluteData(F2, 'vgrid')) < HT_VAR_EPSILON_U), 'Face meshs are not compatible. Not the same vgrid.');
  endif

  if lVMerge
    assert(all(abs(HT_Face_GetAbsoluteData(F1, 'ugrid') - HT_Face_GetAbsoluteData(F2, 'ugrid')) < HT_VAR_EPSILON_U), 'Face meshs are not compatible. Not the same ugrid.');
  endif

  % Intersection was already tested before, but if intersection is not allowed
  % it means there is no intersection and the process can be speed up.
  F1.dims = lDims;
  F1.size = lSize;
  F1.globalPosition = lGlobalPosition;
  F1.nodes = [F1.nodes; F2.nodes];
  F1.pos = lPos;
  assert((isnumeric(F1.r) || isempty(F1.r)) && (isnumeric(F2.r) || isempty(F2.r)), 'Merging to faces with non numeric contact resistance is not supported yet');

  if isempty(F1.r) && isempty(F2.r)
  else
    if isempty(F1.r), F1.r = zeros(numel(F1.nodes), 1); endif;
    if isempty(F2.r), F2.r = zeros(numel(F2.nodes), 1); endif;
    if isscalar(F1.r), F1.r = repmat(F1.r, numel(F1.nodes), 1); endif;
    if isscalar(F2.r), F2.r = repmat(F2.r, numel(F2.nodes), 1); endif;
    F1.r = [F1.r ; F2.r];
  endif

  if ~lAllowIntersect
    if lUMerge, F1.n(1) += F2.n(1); endif;
    if lVMerge, F1.n(2) += F2.n(2); endif;
  else
    [F1.dims, ind] = unique(F1.dims, 'rows');
    F1.nodes = F1.nodes(ind);
    F1.pos = F1.pos(ind);
    if ~isempty(F1.r), F1.r = F1.r(ind); endif;
  endif

endfunction
