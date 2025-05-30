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
% Build one or several node objects.
% To create several nodes at the same time, parameters must be specified using
% cell array mechanism. All parameter does not have to be a cell array, but one
% parameter with a cell array as value defines the number of nodes to be created.
%
% Input arguments: HT_Node_Init(name, parameter, value, ...)
% -> name [string or cell array of strings]
%    First argument is the node name or a cell array of names if several nodes
%    have to be created.
% -> pairs (parameter, value) a set of parameters used to create the nodes
%    -> globalPosition [matrix dim 3x1]: the position of the node in the general
%           axis system. Could be omitted.
%    -> resistance [scalar]: an optinal thermal resistance that will be introduced
%           in each thermal connection with this node
%    -> mode [string]: ['distributed'] or 'single'. Defines how is connected this
%           node to other nodes that will be connected to it.
%           In 'distributed' mode, each thermal connection are independant from each
%           other and node.resistance is introduced in all thermal connection.
%           This mode is used typically with external node "ExternalAirTemperature".
%           In 'single' mode, each thermal connection will share the same thermal
%           resistance there is a slight more complexe connection mechanism. One has
%           indeed to know all objects that will be connected.
%           'copy' [node object] is used to duplicate a node object. All parameters
%           except the name are copied from this object.
%
% Returns:
% A array of structures (of type node)
% ========================================================================
function N = HT_Node_Init(varargin)
  assert(numel(varargin) >= 1, 'At least one parameter is required to init node');
  name = varargin{1};
  
  prop = varargin(2:2:end);
  value = varargin(3:2:end);
  
  % Default parameter
  globalPosition = zeros(3,1);
  resistance = 0;
  mode = 'distributed'; % Distributed means any connection to this node is made parallel to existing connection to this node
                        % In single mode, if a resistance is defined for this node, the whole set of connections must be taken into account.
  
  n = numel(prop);
  for i=1:n
    if strcmpi(prop{i}, "globalPosition")
      globalPosition = value{i};
    elseif strcmpi(prop{i}, "resistance")
      resistance = value{i};
    elseif strcmpi(prop{i}, "mode")
      assert(any(strcmpi(value{i}, {'distributed', 'single'})), 'Invalid value for field <mode>');
      mode = value{i};
    elseif strcmpi(prop{i}, "copy")
      assert(HT_CheckType(value{i}, 'node'), 'Invalid node object');
      obj = value{i};
      globalPosition = obj.globalPosition;
      resistance = obj.resistance;
      mode = obj.mode;
      clear obj;
    else
      error(sprintf('Invalid parameter <%s>', disp(prop{i})));
    endif
  endfor
  
  N = struct( '__type__', 'node', ...
              'name', name, ...
              'globalPosition', globalPosition, ...
              'resistance', resistance, ...
              'mode', mode);
                                         
  for i=1:numel(N)
    assert(numel(N(i).globalPosition) == 3, sprintf('Invalid field <globalPosition> for node <%s>', name));
    assert(numel(N(i).resistance) == 1, sprintf('Invalid field <resistance> for node <%s>', name));
    assert(strcmpi(N(i).mode, 'distributed'), 'Only distributed mode is implemented yet');
    
    N(i).globalPosition = N(i).globalPosition(:);
  endfor
 
endfunction
