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

%======================================================================
%> @brief Returns the intersection of two faces
%>        Warning: the first face has a special effect.
%>        "face1 intersect face 2" is not equal to "face 2 intersect face 1"
%>
%> @param face1 an object of type <face>
%> @param face2 an object of type <face>
%>
%> @retval F a new face
%======================================================================
function F = HT_Face_Intersect(name, Fs, varargin)
  if nargin < 3, varargin = {}; endif;

  assert(ischar(name), 'Face_Intersect::First argument must be a name');
  assert(nargin >= 2, 'Face_Intersect::There must be at least 2 input arguments');
  assert(mod(numel(varargin),2) == 0, 'Face_Intersect::Invalid input parameters');

  prop = varargin(1:2:end);
  values = varargin(2:2:end);

  % Default parameter values
  lParams = struct( 'intersectMesh', false, ...
                    'mustBeCoplanar', false, ...
                    'allowNodeClip', false);
##  lIntersectMesh = false;
##  lMustBeCoplanar = false;
##  lAllowNodeClip = false;

  for i=1:numel(prop)
    if strcmpi(prop{i}, "mesh")
      assert(islogical(values{i}), 'Face_Intersect::Invalid parameter <mesh>');
      lParams.intersectMesh = values{i};
    elseif strcmpi(prop{i}, "coplanar")
      assert(islogical(values{i}), 'Face_Intersect::Invalid parameter <coplanar>');
      lParams.mustBeCoplanar = values{i};
    elseif strcmpi(prop{i}, "nodeclip")
      assert(islogical(values{i}), 'Face_Intersect::Invalid parameter <nodeclip>');
      lParams.allowNodeClip = values{i};
    else
      error(sprintf('Face_Intersect::Invalid input parameter <%s>', prop{i}));
    endif
  endfor

  if isempty(Fs)
    F = {};
    return;
  endif

  if ~isstruct(Fs) && ~iscell(Fs)
    error('Face_Intersect::Invalid array of faces');
  endif

  % Check input faces
  assert(all(arrayfun(@(v) HT_CheckType(v, 'face'), Fs)), 'Invalid input faces');

  % Check orientation
  F = Int_GetObject(Fs(1)); % GetObject returns the object stored in a cell (1x1)
  t = arrayfun(@(v) norm(cross(F.norm, Int_GetObject(v).norm), Inf) < 1E-10, Fs);
  assert(all(t), sprintf('Face_Intersect::Faces <%s> are not oriented in the same direction', strjoin(HT_GetStructField(Fs(t), 'name', true), '; ')));

  % Check coplanarity
  if lParams.mustBeCoplanar
    global HT_VAR_EPSILON_POS

    lPlanePos = dot(F.globalPosition, F.norm);
    t = arrayfun(@(v) abs(dot(Int_GetObject(v).globalPosition, F.norm) - lPlanePos) < HT_VAR_EPSILON_POS, Fs(2:end));
    assert(all(t), ...
            sprintf('Face_Intersect::Faces are not coplanar <%s>', strjoin(HT_GetStructField(Fs(t), 'name', true), '; ')));
    clear lPlanePos;
  endif

  % Lookup all faces and extract the type (general face is .vertex is defined, rectangle based face otherwise)
  lFaceTypeList = arrayfun(@(v) isempty(Int_GetObject(v).vertex), Fs);
  lFacePolygonType = arrayfun(@(v) columns(Int_GetObject(v).dims), Fs);
  lFaceTypeRectangleBased = all(lFaceTypeList);
  lFaceTypeGeneral = all(~lFaceTypeList);

  if lFaceTypeRectangleBased % All faces are rectangle based ?
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
    lIntersection = HT_Rect2.fromIntersection(lRectList);
    assert(~lIntersection.isNull(), 'Face_Intersect::Could not build the union of faces');

    if ~lParams.intersectMesh % only merge the geometry
      F = HT_Face_Init(name,  'size', lIntersection.size(), ...
                              'axis', F.axis, ...
                              'norm', F.norm, ...
                              'globalPosition', lIntersection.to3dSpace(lIntersection.origin(), [F.axis F.norm], F.globalPosition));
    else
      for i=2:numel(Fs)
        F = Int_GetIntersection(F, Int_GetObject(Fs(i)), lParams.allowNodeClip);
      endfor
    endif

  elseif lFaceTypeGeneral % General type face
    lFacePolygonType = unique(lFacePolygonType)
    assert(isscalar(lFacePolygonType), 'Not implemented. Face with different polygon type can not be mixed.');
    assert(lFacePolygonType == 4, 'Not implemented. Only quad based face is supported yet.');
    assert(~lParams.allowNodeClip, 'Not implemented. Only strict match comparison is supported');
    assert(numel(Fs) == 2, 'Not implemented. Only two faces a supported yet.');
    % It is rather a GetSubFace instead of a intersection function.
    % To be changed

    F = Int_IntersectQuadBasedFace(Fs, lParams);
  else
    error('Not implemented. Face type <general/rectangle> can not be mixed');
  endif

endfunction

% Return an object. If obj is a structure, it does nothing.
% If a cell, it returns the content.
function obj = Int_GetObject(obj)
  if iscell(obj)
    obj = obj{1};
  endif
endfunction

function F1 = Int_GetIntersection(F1, F2, lAllowNodeClip)
  global HT_VAR_EPSILON_POS
  global HT_VAR_EPSILON_U

  error('Face_Intersect::Mesh intersection is not yet implemented');
endfunction

function F = Int_IntersectQuadBasedFace(Fs, lParams)
  assert(numel(Fs) == 2);

  F = Int_GetObject(Fs(1));
  F2 = Int_GetObject(Fs(2));

  if lParams.intersectMesh
    % Build the set of Quad2 object for face F
    lFVertex = HT_Face_GetAbsoluteData(F, 'vertex');
    lF2Vertex = HT_Face_GetAbsoluteData(F2, 'vertex');

    lFQuadList = HT_Quad2.from3dSpace(  lFVertex(:, F.dims'(:)), F.axis  );
    lF2QuadList = HT_Quad2.from3dSpace(  lF2Vertex(:, F2.dims'(:)), F.axis  );
    lRemoveQuad = true(rows(F2.dims), 1);

    for i=1:numel(lFQuadList)
      R = lFQuadList{i}.containsQuads(lF2QuadList);
      lRemoveQuad(R) = false;
    endfor

    F2.dims(lRemoveQuad, :) = [];
    F2.nodes(lRemoveQuad) = [];
    F2.pos(lRemoveQuad, :) = [];
    F2.r(lRemoveQuad) = [];
    F = F2;

    % TODO: Polygons and edges are not correctly intersected
    % Now

##    lQuadList = arrayfun(@(v) HT_Quad2.from3dSpace(  lFVertex(F.dims(v,:)), F.axis  ), 1:rows(F.dims), 'uniformOutput', false);

##    lQuadList = arrayfun(@(v) HT_Quad2.from3dSpace( Int_GetObject(v).globalPosition, ...
##                                                    Int_GetObject(v).axis,
##                                                    Int_GetObject(v).size,
##                                                    [F.axis F.norm]), ...
##                         Fs, 'UniformOutput', false);


  else % Intersect only with polygons (larger than meshes)
    error('Not implemented');
  endif
endfunction




