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
%
% This function draws a face or a list of face in a 3D axis
%
% Input arguments are: faces [parameter, value, ...]
% The first input argument is a face or a list of faces.
%
% Next arguments are options, which include:
% .temperature = [vector] a column vector containing temperature value used to
%                compute the appropriate color
% .nodes = [cell of string] containing the name of each node corresponding to
%          the vector of temperature <.temperature>. It is used to retrieve
%          the appropriate temperature component based on node names defined in
%          each face to draw.
% .alpha = [value] transparency coefficient used to draw faces
% .explode = [true/false] Draw faces using an "explosion" mode. It allows different
%          objects to be seen properly. Please note that objects drawn in a single
%          call to Plot_Face function are kept together. An explosition effect
%          occurs only if several calls to Plot_Face function are made.
% .explodeCenter = [vector dim 3x1] Defines the center of the "explosition" operation
%          By default [0,0,0]
% .explodeOffset = [vector dim 3x1] Add an offset to the explosition effect. This
%          may be necessary when small volumes are close to each other.
% .displaygrid = [true/false] If true, the grid, i.e. spatial discretization is
%          displayed for each face. If false, only the face edges are drawn.
% .displayedge = [true/false] Display the edges of each faces if true.
% .edgecolor = [color] the color used to draw edges (ex: 'b' or [0 0 1])
% .facecolortype = ['none'/'flat'/'interp'] select the algorithm used to compute
%                  the face color
% .displayNormal = [true/false] draw the face normal vector
% .displayAxis = [true/false] draw the face normal vector
function hvec = HT_Plot_Face(varargin)
  HT_ImportConstants();

  assert(numel(varargin) >= 1);

  faces = varargin{1};
  if isstruct(faces)
    lOldFaces = faces;
    faces = cell(numel(lOldFaces), 1);
    for k=1:numel(lOldFaces)
      faces{k,1} = lOldFaces(k);
    endfor
    clear lOldFaces;
  endif;

  assert(iscell(faces));
  assert(all(cellfun(@(v) HT_CheckType(v, 'face'), faces)), 'First argument must be a cell array of faces');

  params = struct();
  hvec = [];

  % If the number of input argument is even, it means there
  % is an additionnal structure options
  if mod(numel(varargin), 2) == 0,
    params = varargin{end};
    assert(isstruct(params), "Last argument must be a structure");
  endif

  prop = varargin(2:2:(end-1));
  value = varargin(3:2:end);
  assert(numel(prop) == numel(value), 'Invalid argument list');

  for i=1:numel(prop)
    assert(any(strcmpi(prop{i}, ...
      { 'temperature', 'nodes', 'alpha', 'explode', 'explodeCenter', 'explodeOffset', ...
        'displaygrid', 'displayedge', 'edgecolor', 'gridcolor', 'facecolor', 'facecolortype', 'edgewidth', 'gridwidth', ...
        'displayNormal', 'displayAxis' })), ...
        sprintf('Invalid parameter <%s>', prop{i}));
    params = setfield(params, prop{i}, value{i});
  endfor

  params = HT_CheckField(params, 'temperature', [], {@(v) isvector(v) || isempty(v)});
  params = HT_CheckField(params, 'nodes', [], {@(v) isvector(v) || isempty(v)});
  assert(size(params.temperature, 1) == numel(params.nodes));

  params = HT_CheckField(params, 'alpha', 1.0, {@(v) (v>=0.0) && (v<=1.0)});
  params = HT_CheckField(params, 'explode', 0, {@(v) ((numel(v) == 1) || (numel(v) == 3)) && all(v >= 0)});
  params = HT_CheckField(params, 'explodeCenter', zeros(3,1), {@(v) numel(v) == 3});
  params = HT_CheckField(params, 'explodeOffset', zeros(3,1), {@(v) numel(v) == 3});
  params = HT_CheckField(params, 'displaygrid', true, {@(v) islogical(v) && isscalar(v)});
  params = HT_CheckField(params, 'displayedge', false, {@(v) islogical(v) && isscalar(v)});
  params = HT_CheckField(params, 'edgecolor', 'material', {@(v) (numel(v) == 3) || ischar(v)});
  params = HT_CheckField(params, 'gridcolor', 'material', {@(v) (numel(v) == 3) || ischar(v)});
  params = HT_CheckField(params, 'facecolor', 'material', {@(v) (numel(v) == 3) || ischar(v)});
  params = HT_CheckField(params, 'facecolortype', 'none', {@(v) any(strcmpi(v, {'none','flat','interp'}))});
  params = HT_CheckField(params, 'edgewidth', 1, {@(v) numel(v) == 1});
  params = HT_CheckField(params, 'gridwidth', 1, {@(v) numel(v) == 1});
  params = HT_CheckField(params, 'displayNormal', false, {@(v) islogical(v) && isscalar(v)});
  params = HT_CheckField(params, 'displayAxis', false, {@(v) islogical(v) && isscalar(v)});

  if numel(params.explode) == 1,
    params.explode = repmat(params.explode, 3, 1);
  endif

  explodeVec = zeros(3,1);

  % Apply explode option
  if any(params.explode > 0)
    for i=1:numel(faces)
      explodeVec += faces{i}.globalPosition + 0.5 * faces{i}.axis * faces{i}.size(:);
    endfor
    explodeVec = (explodeVec/max(numel(faces), 1) - params.explodeCenter(:)) .* params.explode(:);
    explodeVec += params.explodeOffset(:);
  endif

  if params.displaygrid
    for i=1:numel(faces)
      face = faces{i};
      nf = numel(face.nodes);

      % If no mesh, no nodes is present, continue
      if nf == 0
        continue;
      endif

      dims = face.dims;

      if ~isempty(params.temperature) && ~isempty(params.nodes)
        C = HT_Result_GetNodeTemperature(face.nodes, params.temperature, params.nodes);
        dims(isnan(C),:) = [];
        C(isnan(C)) = [];
        nf = rows(dims);
      elseif ~ischar(params.facecolor)
        C = repmat(params.facecolor(:)', nf, 1);
      elseif strcmpi(params.facecolor, 'material') && ~isempty(face.material)
        if numel(face.material) == 1
          C = face.material.color(:)';
##          C = repmat(face.material.color(:)', nf, 1);
        else
          assert(~isempty(face.materialIndex), 'Invalid field <materialIndex>');
          C = cell2mat(arrayfun(@(mat) mat.color(:)', face.material(:), 'Uniformoutput', false));
          C = C(face.materialIndex,:);
        endif
      else
        C = ones(nf, 3);
      endif

      % If .vertex is not defined, it means nodes are squares
      if isempty(face.vertex)
        dims = [       (dims(:,[1 1 2 2])')(:)' ; ...
                       (dims(:,[3 4 4 3])')(:)'];

        dims .*= face.size(:); % Multiply each row by the size

        % TODO: remove redundant vertices
        V = face.globalPosition + face.axis * dims;

        VFace = reshape(1:(4*nf), 4, nf)';
      else
        V = HT_Face_GetAbsoluteData(face, 'vertex');
        VFace = face.dims;
      endif

      V = round(V'/HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS; % Arrondis au 1/10mm

      lGridColor = params.gridcolor;

      if strcmpi(params.gridcolor, 'material') && (numel(face.material) == 1)
        lGridColor = 0.5*face.material.color;
      else
        lGridColor = [0, 0, 0];
      endif

      V += explodeVec';

      h = patch("Faces", VFace, "Vertices", V, ...
##                            "CData", repmat(C', 4, nf), ...
                            "FaceVertexCData", C, ...
                            "facecolor", params.facecolortype, ...
                            "facealpha", params.alpha,...
                            "edgecolor", lGridColor, ... 'none', ...lGridColor, ...
                            "linewidth", params.gridwidth, ...
                            "facelighting", "none");
      hvec = [hvec, h];
    endfor
  endif

  if params.displayedge
    for i=1:numel(faces)
      face = faces{i};

      lEdgeColor = params.edgecolor;
      lFaceColor = params.facecolor;

      if strcmpi(params.edgecolor, 'material') && (numel(face.material) == 1)
        lEdgeColor = face.material.color;
      elseif ~isnumeric(params.edgecolor) || (numel(params.edgecolor) != 3)
        lEdgeColor = [0, 0, 0];
      endif

      if params.displaygrid
        lFaceColor = [];
      elseif strcmpi(lFaceColor, 'material') && ~isempty(face.material)
        if numel(face.material) == 1
          C = face.material.color(:)';
        else
          assert(~isempty(face.materialIndex), 'Invalid field <materialIndex>');
          C = cell2mat(arrayfun(@(mat) mat.color(:)', face.material(:), 'Uniformoutput', false));
          C = C(face.materialIndex,:)
        endif
      elseif ~isnumeric(lFaceColor) || (numel(lFaceColor) != 3)
        lFaceColor = [1, 1, 1];
      endif

      if ~isempty(face.polygons)
        V = HT_Face_GetAbsoluteData(face, 'vertex');
        V += explodeVec';
        V = round(V/HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS; % Arrondis au 1/10mm

        if isempty(lFaceColor)
##          for k=1:rows(face.edges)
##            h = line(   "xdata", V(1,face.edges(k,:)), ...
##                        "ydata", V(2,face.edges(k,:)), ...
##                        "zdata", V(3,face.edges(k,:)), ...
##                        "linewidth", params.edgewidth, ...
##                        "color", lEdgeColor);
##          endfor

        else
          h = patch("Faces", face.polygons, "Vertices", V', ...
                                    "FaceVertexCData", lFaceColor, ...
                                    "FaceColor", params.facecolortype, ...
                                    "linewidth", params.edgewidth, ...
                                    "edgecolor", lEdgeColor, ...
                                    "facelighting", "none");
        endif
      else % Old version with quads only
        if isempty(lFaceColor)
          C = HT_Face_GetAbsoluteData(face, 'corners');
          h = line(   "xdata", [C(1,:), C(1,1)], ...
                      "ydata", [C(2,:), C(2,1)], ...
                      "zdata", [C(3,:), C(3,1)], ...
                      "linewidth", params.edgewidth, ...
                      "color", lEdgeColor);
        else
          C = HT_Face_GetAbsoluteData(face, 'corners');
          h = patch("Faces", [1 2 3 4], "Vertices", C', ...
                                    "FaceVertexCData", lFaceColor, ...
                                    "FaceColor", params.facecolortype, ...
                                    "linewidth", params.edgewidth, ...
                                    "edgecolor", lEdgeColor, ...
                                    "facelighting", "none");
        endif
      endif

      hvec = [hvec, h];
    endfor
  endif

  if params.displayNormal
    for i=1:numel(faces)
      face = faces{i};
      C = HT_Face_GetAbsoluteData(face, 'center');
      A = face.norm * mean(face.size);
##      drawArrow3d(C, face.norm, [0 1 0]);
      h = quiver3(C(1), C(2), C(3), A(1), A(2), A(3));
      set(h, 'linewidth', 2.0, 'color', [1 0 1]);
    endfor
  endif

  if params.displayAxis
    for i=1:numel(faces)
      face = faces{i};
      C = face.globalPosition;
      A = face.axis(:,1) * mean(face.size);
##      drawArrow3d(C, face.norm, [0 1 0]);
      h = quiver3(C(1), C(2), C(3), A(1), A(2), A(3));
      set(h, 'linewidth', 2.0, 'color', [1 0 0]);
      A = face.axis(:,2) * mean(face.size);
      h = quiver3(C(1), C(2), C(3), A(1), A(2), A(3));
      set(h, 'linewidth', 2.0, 'color', [0 1 0]);
    endfor
  endif
endfunction

