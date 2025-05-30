%  This file is part of project HiveTemp.
%  This work was supported by the Better-B project, which has received funding
%  from the European Union, the Swiss State Secretariat for Education, Research
%  and Innovation (SERI) and UK Research and Innovation (UKRI) under the UK
%  government's Horizon Europe funding guarantee (grant number 10068544). Views
%  and opinions expressed are however those of the author(s) only and do not
%  necessarily reflect those of the European Union, European Research Executive
%  Agency (REA), SERI or UKRI. Neither the European Union nor the granting
%  authorities can be held responsible for them.
%
%  Copyright (c) 2023-2025 AltRD-Emmanuel Ruffio
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
% This function returns an empty model structure
%
% Input arguments:
% <type> = [string] or {cell of string} specifying the model type and subtype
% <name> = [string] the name of the model
% <params> = [struct] data stucture ('paramname', paramvalue, ...)
%            containing parameters
%
function M = HT_Model_Hygro_Init(type, name, params)
  assert(nargin >= 1, 'model type must be defined');

  if nargin < 2, name = ''; endif;
  if nargin < 3, params = struct(); endif;

  if (iscell(type)) % Add the type 'model'
    if size(type,2) > 0
      type = unique([type, 'model_hygro']);
    else
      type = unique([type; 'model_hygro']);
    endif
  elseif (ischar(type))
    type = {'model_hygro', type};
  else
    error('Argument <type> is not valid');
  endif

##  M.nodes = [lLayerData.nodes](:);   % Extract and merge all nodes from all layers
##  M.materials = [lLayerData.material](:); % Material for each layer
##  M.nodesMat = vertcat(lLayerData.nodesMat);    % Extract and merge all conductance matrix based on geometry only
##  M.init0Vec = [lLayerData.init0];       % Initialization value
##  M.volumeVec = [lLayerData.volumeVec]; % Volume of each nodes
##  M.globalPosition = lParams.globalPosition;
##  M.axis = lParams.axis;

  M = struct( 'name', name, ...
              '__type__', [],...
              'CVec', [],          ...   % [vector] Capacity of each nodes
              'nodes', [], ...
              'nodesMat', [], ...     % [dim Nx13] Relation between nodes
              'materials', [], ...
              'init0Vec', [], ...     % [dim Nnodes x 1] Initial value
              'volumeVec', [], ...     % [dim Nnodes x 1] Volume of each nodes
              'rad', {{}}, ...          % Radiation models
              'params', params, ...   % Store the parameters used to create the model
              'globalPosition', [], ...
              'axis', [], ...
              'submodel', [], ...
              'timer', [], ...
              'counter', []);

  % .nodesMat contains:
  % IndexNode1 IndexNode2 IndexMat1 IndexMat2 Vec1(node to border) Vec2(border to node) RGeometry1 RGeometry2 Rcontact

  M.__type__ = type;  % Type may be a cell, since a structure may have severals types
##  M.rad = {}; % Can not initialize with an empty cell since struct interpet this as a struct array
  M.timer = struct('buildingTime', 0, 'solvingTime', 0);
  M.counter = struct('warnings', 0, 'errors', 0);
endfunction
