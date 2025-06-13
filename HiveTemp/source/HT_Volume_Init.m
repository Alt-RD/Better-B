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

function F = HT_Volume_Init(varargin)
  global HT_VAR_EPSILON_POS;
  assert(numel(varargin) >= 1, 'At least one parameter is required to init face');

  lStructMode = isstruct(varargin{1}); % iscell(varargin{1}) ||;

  % Default parameter
  name = {};
  shapeType = 'undefined';
  dimension = [];
  nodes = [];
  nodesVolume = [];

  if ~lStructMode
    prop = varargin(1:2:end);
    value = varargin(2:2:end);

    assert(numel(prop) == numel(value), 'Invalid input parameters');

##    round([lPos1; lPos2] ./ HT_VAR_EPSILON_POS) * HT_VAR_EPSILON_POS

    nProp = numel(prop);
    for i=1:nProp
      if strcmpi(prop{i}, "name")
        name = value{i}; %HT_Round(value{i}, HT_VAR_EPSILON_POS);
      elseif strcmpi(prop{i}, "shapeType")
        shapeType = value{i}; %HT_Round(value{i}, HT_VAR_EPSILON_POS);
      elseif strcmpi(prop{i}, "dimension")
        dimension = value{i};
      elseif strcmpi(prop{i}, "nodes")
        nodes = value{i};
      elseif strcmpi(prop{i}, "nodesVolume")
        nodesVolume = value{i};
      elseif strcmpi(prop{i}, "copy")
        assert(HT_CheckType(value{i}, 'volume'), 'Invalid face object');
        f = value{i};
        name = f.name;
        shapeType = f.shapeType;
        dimension = f.dimension;
        nodes = f.nodes;
        nodesVolume = f.nodesVolume;
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

  F = struct( '__type__', 'volume', ...
              'name', name(:), ...
              'shapeType', shapeType(:), ...
              'dimension', dimension(:), ...
              'nodes', nodes(:), ...
              'nodesVolume', nodesVolume(:));
    % Additionnal resistance to add for face connection
    % May be a single scalar (Km²/W) or a vector defining the resistance for each nodes

  for i=1:numel(F)
    F(i).dimension = F(i).dimension(:);
    F(i).nodes = F(i).nodes(:);
    F(i).nodesVolume = F(i).nodesVolume(:);

##    assert(isempty(F(i).materialIndex) || ((numel(F(i).materialIndex) == numel(F(i).nodes)) && all((F(i).materialIndex > 1) | (F(i).materialIndex <= numel(F(i).material)))), 'Invalid parameter <materialIndex>');
##    assert(ischar(F(i).model), 'Invalid parameter <model>');
  endfor

endfunction
