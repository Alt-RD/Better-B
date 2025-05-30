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
%> @brief Brief Create a circular mesh from an existing face
%>        Clear the existing mesh if necessary
%>
%> @param face an object of type <face>
%> Property:
%> -> rGrid: [vector of double] in [0;1]: boundaries of nodes. Ex: [0, 0.1, 0.2, 1]
%>           [integer]: create a uniformly spaced nodes along the radius
%>
%> @retval F a new face
%======================================================================
function Fs = HT_Face_CreateCircularMesh(Fs, varargin)
  global HT_VAR_EPSILON_POS;
  global HT_VAR_EPSILON_U;

  assert(nargin >= 1, "Missing argument");
  assert(mod(numel(varargin), 2) == 0, "Invalid input parameters");

  lParams = struct( 'radius', [], ...
                    'theta', [0, 2*pi], ...
                    'rGrid', [0, 1], ...
                    'thetaGrid', int32(8), ...
                    'globalPosition', NA(3,1), ...
                    'size', NA(2,1), ...
                    'axis', NA(3,2), ...
                    'normal', NA(3,1));

  % If Fs is a face based on circular type, parameters are extracted first from
  % this face
  if HT_CheckType(Fs, 'face') && ~isempty(Fs.metadata) && strcmpi(Fs.metadata.type, 'circle')
    lParams.radius = Fs.metadata.uSize;
    lParams.theta = Fs.metadata.vSize;
    lParams.rGrid = Fs.metadata.uGrid;
    lParams.thetaGrid = Fs.metadata.vGrid;
  endif

  prop = varargin(1:2:end);
  value = varargin(2:2:end);

  for i=1:numel(prop)
    if strcmpi(prop{i}, 'radius')
      lParams.radius = value{i};
    elseif strcmpi(prop{i}, 'theta')
      lParams.theta = value{i};
    elseif strcmpi(prop{i}, 'rGrid')
      lParams.rGrid = value{i};
    elseif strcmpi(prop{i}, 'thetaGrid')
      lParams.thetaGrid = value{i};
    elseif strcmpi(prop{i}, 'globalPosition')
      lParams.globalPosition = value{i};
    elseif strcmpi(prop{i}, 'size')
      lParams.size = value{i};
    elseif strcmpi(prop{i}, 'axis')
      lParams.axis = value{i};
    elseif strcmpi(prop{i}, 'normal')
      lParams.normal = value{i};
    else
      error(sprintf('Invalid input parameter name <%s>', prop{i}));
    endif
  endfor

  if isempty(Fs)
    Fs = HT_Face_Init('null');
  endif

  Fs.dims = [];
  Fs.vertex = [];
  Fs.edges = [];
  Fs.polygons = [];

  if ~isempty(lParams.size),            Fs.size = HT_MergeArray(Fs.size, lParams.size, 2); endif
  if ~isempty(lParams.normal),          Fs.norm = HT_MergeArray(Fs.norm, lParams.normal, 2); endif
  if ~isempty(lParams.axis),            Fs.axis = HT_MergeArray(Fs.axis, lParams.axis, 2); endif
  if ~isempty(lParams.globalPosition),  Fs.globalPosition = HT_MergeArray(Fs.globalPosition, lParams.globalPosition, 2); endif

  assert(~any(isna(Fs.size)), 'Parameter <size> is not valid');
  assert(~any(isna(Fs.norm)), 'Parameter <normal> is not valid');
  assert(~any(isna(Fs.axis)), 'Parameter <axis> is not valid');
  assert(~any(isna(Fs.globalPosition)), 'Parameter <globalPosition> is not valid');

  if isempty(lParams.radius), lParams.radius = 0.5*Fs.size(1); endif

  assert(~any(isna(lParams.radius)) && ~isempty(lParams.radius), 'Parameter <radius> is not defined');

  % Create the mesh
  if isinteger(lParams.thetaGrid)
    lAngleVec = double(0:(lParams.thetaGrid-1))'/double(lParams.thetaGrid);
  else
    lAngleVec = lParams.thetaGrid(:);
  endif

  lAngleVec = lParams.theta(1) + lAngleVec * (lParams.theta(2) - lParams.theta(1));

  lVertexUnit = [cos(lAngleVec) sin(lAngleVec)];
  lVertexUnitLoop = [lVertexUnit(end,:) ; lVertexUnit];

  if isinteger(lParams.rGrid),
    lParams.rGrid = double(0:lParams.rGrid)'/double(lParams.rGrid);
  endif

  lParams.rGrid = round(lParams.rGrid / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  assert((lParams.rGrid(1) == 0.0) && (lParams.rGrid(end) == 1.0), 'Invalid <rGrid> parameter');

  nr = numel(lParams.rGrid)-1;
  ny = numel(lAngleVec);

  rDimTotal = lParams.radius(1) + lParams.rGrid * (lParams.radius(end) - lParams.radius(1));

  Fs.vertex = repmat(rDimTotal, ny, 1) .* repelem(lVertexUnit, (nr+1), 1);

  Fs.vertex += 0.5*Fs.size';
  Fs.vertex ./= Fs.size';

  assert(all(Fs.vertex(:,1) >= -1E-10) && all(Fs.vertex(:,1) <= 1+1E-10), 'Discrepancy between <size> and <radius>');
  assert(all(Fs.vertex(:,2) >= -1E-10) && all(Fs.vertex(:,2) <= 1+1E-10), 'Discrepancy between <size> and <radius>');

  Fs.vertex(Fs.vertex < 0) = 0;
  Fs.vertex(Fs.vertex > 1) = 1;


  Fs.dims = int32( (1:(nr*ny))' + repelem((0:(ny-1))',nr,1) + ...
                              [ 0, -nr-1, -nr, 1]); %Vertices are defined the direct order to match the normal vector
  Fs.dims(1:nr,2:3) += (nr+1)*ny;
  Fs.edges = int32([(1:ny)*(nr+1)    , nr+1 ; ... External lines
                           (1:ny)*(nr+1)-nr , 1 ]); ... Internal lines
  Fs.polygons = int32( 1 + (0:(ny-1))'*(nr+1) + ...
                              [ 0, -nr-1, -1, nr]); ... %Vertices are defined the direct order to match the normal vector
  Fs.polygons(1,2:3) += (nr+1)*ny;

##  faces.top.dims = int32( (1:(nr*ny))' + repelem((0:(ny-1))',nr,1) + ...
##                              [ 0, -nr-1, -nr, 1]); %Vertices are defined the direct order to match the normal vector
##  faces.top.dims(1:nr,2:3) += (nr+1)*ny;
##  faces.top.edges = int32([(1:ny)*(nr+1)    , nr+1 ; ... External lines
##                           (1:ny)*(nr+1)-nr , 1 ]); ... Internal lines
##  faces.top.polygons = int32( 1 + (0:(ny-1))'*(nr+1) + ...
##                              [ 0, -nr-1, -1, nr]); ... %Vertices are defined the direct order to match the normal vector
##  faces.top.polygons(1,2:3) += (nr+1)*ny;

endfunction


