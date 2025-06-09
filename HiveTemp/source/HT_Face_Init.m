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
% Build a or several face objects
%
% Input arguments: pairs of 'parameter name' <-> 'parameter value'
% Parameter value can be a single object (integer, array...) or a cell array.
% If the input parameter is a cell array, several faces will be created and returned.
% 1) 'size' : the dimensions of the face along the (axis(1), axis(2))
% 2) 'nodes': the names of nodes that are related to this face
% 3) 'position': (dim Nx2) relative positions (u,v) of the center of every nodes
% 4) 'dims': (dim Nx4) relative positions/sizes (u1,u2,v1,v2) of every nodes.
%            Ex: [0 1 0 1] of only 1 node is defined for the entire face
%            (dim Nx8) relative positions (u1,v1,u2,v2,u3,v3,u4,v4) for quads
%                      In that cases <.n> is not 1x2 but only 1x1 (number of nodes)
%			 (dim Nx6) for triangles
%                      In that cases <.n> is not 1x2 but only 1x1 (number of nodes)
% 5) 'n': (dim 1x2) number of nodes along U,V direction
% 6) 'axis': (dim 3x2) global vectors corresponding to local axis (u,v)
% 7) 'norm': (dim 3x1) normal vector pointing outside the medium
% 8) 'globalposition': (dim 3x1) position in the global axis system of (u=0, v=0)
% 9) 'r': (dim Nx1) Additionnal resistance to add to each node for face connection
%                  May be a single scalar (Km²/W) or a vector defining the
%                  resistance for each nodes
% 10) 'model': (string) the model name from which the face is derived
% 11) 'copy': (face object) specify a face that will used to retrieve geometry
%             ie. 'size', 'position', 'dims', 'axis', 'norm', 'globalposition'
%             Geometry can still be overwritten by explicitely specifying parameters
%             after the parameter 'base'
%
%
% Output parameters:
% 1) F: a face struct or a face struct array (if some input arguments are cell array)
%
% Changed:
% 15/06/2021: Add normal vector <norm> which points out to the exterior side
%             This information is usefull since local axis u,v may be flipped
% 08/01/2023: Add parameter 'base' to copy parameters from a existing face

%======================================================================
%> @brief Create a new or a copy of a face
%>
%> This function create a new face object from scratch or from a existing face.
%> It supports severals parameters to set or modify properties of the new face.
%> If a parameter is specified several times, the last value overwrites the
%> previous ones. It is useful when one wants to copy a face but still wants
%> some parameters to be changed.
%>
%> @param name Nom de la face
%> @param varargin A list of parameters. Currently supported parameters are:
%>        "size", "nodes", "position", "dims", "n", "axis", "norm",
%>        "globalPosition", "r", "model", "copy".
%>
%>
%> @retval F return value for the first output variable
%======================================================================
function F = HT_Face_Init(varargin)
  global HT_VAR_EPSILON_POS;
  assert(numel(varargin) >= 1, 'At least one parameter is required to init face');

  lStructMode = isstruct(varargin{1}); % iscell(varargin{1}) ||;

  % Default parameter
  faceSize = zeros(2,1);
  nodes = [];
  pos = [];
  dims = int32([]);
  edges = int32([]);    % Matrix of indices that defines the edges of the faces
  polygons = int32([]); % Matrix (Nx3 ou Nx4) that defines the polygons of the solid
  axis = [];
  norm = [];
  n = [];
  globalPosition = zeros(3,1);
  resistance = [];
  vertex = [];
  material = [];
  materialIndex = [];
  model = 'null';
  metadata = struct();

  if ~lStructMode
    name = varargin{1};

    prop = varargin(2:2:end);
    value = varargin(3:2:end);

    assert(numel(prop) == numel(value), 'Invalid input parameters');

##    round([lPos1; lPos2] ./ HT_VAR_EPSILON_POS) * HT_VAR_EPSILON_POS

    nProp = numel(prop);
    for i=1:nProp
      if strcmpi(prop{i}, "size")
        faceSize = HT_Round(value{i}, HT_VAR_EPSILON_POS);
      elseif strcmpi(prop{i}, "nodes")
        nodes = value{i};
      elseif strcmpi(prop{i}, "position")
        pos = value{i};
      elseif strcmpi(prop{i}, "dims")
        dims = value{i};
      elseif strcmpi(prop{i}, "edges")
        edges = value{i};
      elseif strcmpi(prop{i}, "polygons")
        polygons = value{i};
      elseif strcmpi(prop{i}, "n")
        n = value{i};
      elseif strcmpi(prop{i}, "axis")
        axis = value{i};
      elseif strcmpi(prop{i}, "norm")
        norm = value{i};
      elseif strcmpi(prop{i}, "globalPosition")
        globalPosition = HT_Round(value{i}, HT_VAR_EPSILON_POS);
      elseif strcmpi(prop{i}, "r")
        resistance = value{i};
      elseif strcmpi(prop{i}, "model")
        model = value{i};
      elseif strcmpi(prop{i}, "vertex")
        vertex = value{i};
      elseif strcmpi(prop{i}, "material")
        material = value{i};
      elseif strcmpi(prop{i}, "materialIndex")
        materialIndex = value{i};
      elseif strcmpi(prop{i}, "metadata")
        metadata = value{i};
      elseif strcmpi(prop{i}, "copy")
        assert(HT_CheckType(value{i}, 'face'), 'Invalid face object');
        f = value{i};
        nodes = f.nodes;
        faceSize = f.size;
        pos = f.pos;
        dims = f.dims;
        edges = f.edges;
        polygons = f.polygons;
        axis = f.axis;
        norm = f.norm;
        globalPosition = f.globalPosition;
        resistance = f.r;
        vertex = f.vertex;
        model = f.model;
        material = f.material;
        materialIndex = f.materialIndex;
        metadata = f.metadata;
        clear f;
      else
        error(sprintf('Invalid parameter <%s>', disp(prop{i})));
      endif
    endfor
  elseif isstruct(varargin{1})
    lObjectList = varargin{1};
    nObject = numel(lObjectList);

    name = repmat({ '' }, nObject, 1);

    prop = fieldnames(lObjectList);

    nProp = numel(prop);
    for i=1:nProp
      if strcmpi(prop{i}, "name")
        name = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "size")
        faceSize = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "nodes")
        nodes = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "position")
        pos = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "dims")
        dims = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "edges")
        edges = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "polygons")
        polygons = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "n")
        n = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "axis")
        axis = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "norm")
        norm = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "globalPosition")
        globalPosition = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "r")
        resistance = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "vertex")
        vertex = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "model")
        model = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "material")
        material = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "materialIndex")
        materialIndex = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "metadata")
        metadata = arrayfun(@(v) v.(prop{i}), lObjectList, 'UniformOutput', false);
      elseif strcmpi(prop{i}, "copy")
        assert(numel(v.(prop{i})) == 1, 'Invalid face object. Must be a single object. Not an array');
        assert(HT_CheckType(v.(prop{i}), 'face'), 'Invalid face object');
        nodes = arrayfun(@(v) v.copy.nodes, lObjectList, 'UniformOutput', false);
        faceSize = arrayfun(@(v) v.copy.faceSize, lObjectList, 'UniformOutput', false);
        pos = arrayfun(@(v) v.copy.pos, lObjectList, 'UniformOutput', false);
        dims = arrayfun(@(v) v.copy.dims, lObjectList, 'UniformOutput', false);
        edges = arrayfun(@(v) v.copy.edges, lObjectList, 'UniformOutput', false);
        polygons = arrayfun(@(v) v.copy.polygons, lObjectList, 'UniformOutput', false);
        axis = arrayfun(@(v) v.copy.axis, lObjectList, 'UniformOutput', false);
        norm = arrayfun(@(v) v.copy.norm, lObjectList, 'UniformOutput', false);
        globalPosition = arrayfun(@(v) v.copy.globalPosition, lObjectList, 'UniformOutput', false);
        resistance = arrayfun(@(v) v.copy.r, lObjectList, 'UniformOutput', false);
        vertex = arrayfun(@(v) v.copy.vertex, lObjectList, 'UniformOutput', false);
        model = arrayfun(@(v) v.copy.model, lObjectList, 'UniformOutput', false);
        material = arrayfun(@(v) v.copy.material, lObjectList, 'UniformOutput', false);
        materialIndex = arrayfun(@(v) v.copy.materialIndex, lObjectList, 'UniformOutput', false);
        metadata = arrayfun(@(v) v.copy.metadata, lObjectList, 'UniformOutput', false);
        clear f;
      else
        error(sprintf('Invalid parameter <%s>', disp(prop{i})));
      endif
    endfor
  else
    error('Not implemented yet');
  endif

  if ~iscell(nodes) || all(cellfun(@(v) ischar(v), nodes))
    % If all elements of nodes are char, convert it to a single element cell
    % to prevent "struct initialisation" to build Nnodes faces.
    nodes = { nodes };
  endif

  F = struct( '__type__', 'face', ...
              'name', name, ...
              'size', faceSize, ...
              'globalPosition', globalPosition, ...
              'nodes', nodes, ...
              'pos', pos, ...
              'dims', dims, ...
              'edges', edges, ...
              'polygons', polygons, ...
              'n', n, ...
              'axis', axis, ...
              'norm', norm, ...
              'r', resistance, ...
              'material', material, ...
              'materialIndex', materialIndex, ...
              'vertex', vertex,...
              'model', model,...
              'metadata', metadata); % Additionnal resistance to add for face connection
                               % May be a single scalar (Km²/W) or a vector defining the resistance for each nodes

  for i=1:numel(F)
    F(i).size = F(i).size(:);
    F(i).globalPosition = F(i).globalPosition(:);
    F(i).r = F(i).r(:);

    assert(isempty(F(i).nodes) || any(numel(F(i).r) == [0 1 numel(F(i).nodes)]), sprintf('Invalid parameter <r>. Size is %d instead of %d.', numel(F(i).r), numel(F(i).nodes)));
    assert(any(size(F(i).dims, 1) == [0 numel(F(i).nodes)]), 'Invalid parameter <dims>');
    assert(any(size(F(i).dims, 2) == [0 4]), 'Invalid parameter <dims>');
    assert(isempty(F(i).polygons) || (max(max(F(i).polygons)) <= rows(vertex)), 'Invalid parameter <dims>. Some indices exceeds the number of vertices');
    assert(isempty(F(i).polygons) || ~isempty(F(i).vertex), 'Invalid parameter <polygons>. <vertex> are not defined');
    assert(any(numel(F(i).n) == [0 2]), 'Invalid parameter <n>');
    assert(any(size(F(i).pos, 1) == [0 numel(F(i).nodes)]), 'Invalid parameter <pos>');
    assert(any(size(F(i).pos, 2) == [0 2]), 'Invalid parameter <pos>');
    assert(all(size(F(i).axis) == [3 2]) || isempty(F(i).axis), 'Invalid parameter <axis>');
    assert(all(size(F(i).norm) == [3 1]) || isempty(F(i).norm), 'Invalid parameter <axis>');
    assert(all(size(F(i).globalPosition) == [3 1]), 'Invalid parameter <globalPosition>');
    assert(all(size(F(i).size) == [2 1]), 'Invalid parameter <size>');
    assert(isinteger(F(i).dims) || isempty(F(i).vertex), 'Invalid type of <dims>. Must be integer if <vertex> is specified');
    assert(isempty(F(i).vertex) || isempty(F(i).dims) || (max(max(F(i).dims)) <= rows(F(i).vertex)), 'Invalid parameter <dims>. Some indices exceeds the number of vertices');
    assert(isempty(F(i).edges) || ~isempty(F(i).vertex), 'Invalid parameter <edges>. <vertex> are not defined');
    assert(isempty(F(i).material) || HT_CheckType(F(i).material, 'material'), 'Invalid parameter <material>');
    assert(isempty(F(i).materialIndex) || ((numel(F(i).materialIndex) == numel(F(i).nodes)) && all((F(i).materialIndex > 1) | (F(i).materialIndex <= numel(F(i).material)))), 'Invalid parameter <materialIndex>');
    assert(ischar(F(i).model), 'Invalid parameter <model>');
  endfor

endfunction
