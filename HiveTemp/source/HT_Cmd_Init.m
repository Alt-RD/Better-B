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
%  Copyright (c) 2022 Montpellier-University
%  Copyright (c) 2022-2025 AltRD: Emmanuel Ruffio
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
% Build one or several commands.
% Each time the property "name" is found, a new command named <name> is created
% and added to the list of command returned by this function.
%
% This function takes a list of pairs (property, value).
% List of available properties:
% 1) "name" : the command name
% 2) "type" : the command type: temperature or flux
% 3) "data" : the data associated to the command
% 4) "nodes": the list of nodes controled by the command
% 5) "area": the area or any vector of coefficient useful for the command
%
% d
function F = HT_Cmd_Init(varargin)
  assert(numel(varargin) > 0, 'Missing arguments');
  assert(mod(numel(varargin), 2) == 0);

  prop = varargin(1:2:end);
  value = varargin(2:2:end);

  lCmdDefault = struct( '__type__', 'command', ...
                        'name', '', ...
                        'type', '', ...
                        'nodes', {{}}, ...
                        'data', [], ...
                        'area', [], ...
                        'force', false);
  lCurrentCmd = lCmdDefault;
  lCommandNameList = {};
  lCmdList = {};

  for k=1:numel(prop)
    if strcmpi(prop{k}, 'name')
      if ~isempty(lCurrentCmd.name) % Add the previous command object to the list
        lCmdList = [lCmdList; INT_Finalize(lCurrentCmd)];
      endif

      lCurrentCmd = lCmdDefault;

      assert(~isempty(value{k}) && ischar(value{k}), 'Invalid command name');
      assert(~any(cellfun(@(v) strcmpi(value{k}, v), lCommandNameList)), 'Duplicate command name');
      lCurrentCmd.name = value{k};
      lCommandNameList = [lCommandNameList, value{k}];

    elseif strcmpi(prop{k}, 'nodes')
      assert( ~isempty(value{k}), sprintf('Invalid node name/object'));

      if iscell(value{k}) && all(cellfun(@(v) ischar(v), value{k}))
        assert(all(cellfun(@(v) ischar(v) && ~isempty(v), value{k})), 'Invalid field <nodes>');
        lCurrentCmd.nodes = value{k};
      elseif ischar(value{k})
        lCurrentCmd.nodes = { value{k} };
      elseif all(HT_CheckType(value{k}, 'face'))
        lCurrentCmd.nodes = value{k};
      else
        error('Invalid field <nodes>');
      endif

    elseif strcmpi(prop{k}, 'type')
      assert(any(strcmpi(value{k}, {'temperature', 'flux'})), 'Invalid command field <type>');
      lCurrentCmd.type = value{k};

    elseif strcmpi(prop{k}, 'data')
      if ischar(value{k})
        assert(exist(value{k}, 2), sprintf('Invalid file path <%s>', value{k}));
      elseif is_function_handle(value{k})
        error('Not implemented');
      elseif isreal(value{k})
        lCurrentCmd.data = value{k};
      elseif iscell(value{k})
        lCurrentCmd.data = value{k};
      else
        error('Invalid field <data>');
      endif

    elseif strcmpi(prop{k}, 'area')
      lCurrentCmd.area = value{k};
    elseif strcmpi(prop{k}, 'force')
      lCurrentCmd.force = value{k};
    else
      error(sprintf('Invalid parameter <%s>', strtrim(disp(prop{k}))));
    endif
  endfor

  if ~isempty(lCurrentCmd.name) % Add the previous command object to the list
    lCmdList = [lCmdList; INT_Finalize(lCurrentCmd)];
  endif

  F = lCmdList;

endfunction

function C = INT_Finalize(C)
  HT_ImportConstants();

  if any(HT_CheckType(C.nodes, 'face'))
    assert(all(HT_CheckType(C.nodes, 'face')), 'Invalid <nodes>. Must be objects of type <face>');

    if isstruct(C.nodes) % Convert struct array to cell array
      C.nodes = arrayfun(@(v) v, C.nodes, 'UniformOutput', false);
    endif

    assert(iscell(C.nodes), 'Invalid parameter <nodes>');

    % Check <data>
    if ~iscell(C.data) % One cell for each face object
      assert(~isempty(C.data), 'No data specified');
##      assert(numel(C.nodes) == 1, 'A cell array must be provided for <data> when using a cell array of faces');
      C.data = repmat({ C.data }, numel(C.nodes), 1);
    endif

    assert(numel(C.data) == numel(C.nodes), 'Invalid size of parameter <data>');
    t_data = cellfun(@(v) rows(v), C.data);
    t_nodes = cellfun(@(v) numel(v.nodes), C.nodes);
    assert(all(arrayfun(@(i) any(t_data(i) == [1, t_nodes(i)]), 1:numel(t_data))), 'The size of some data matrix does not match the number of nodes');

    % Check <area>
    if ischar(C.area)
      assert(strcmpi(C.area, 'none'), 'Invalid field <area>');
      C.area = repmat({ [] }, numel(C.nodes), 1);
    elseif isempty(C.area)
      C.area = HT_Face_GetAbsoluteData(C.nodes, 'nodesArea');
    elseif ~iscell(C.area)
##      assert(numel(C.nodes) == 1, 'A cell array must be provided for <data> when using a cell array of faces');
      C.area = repmat({ C.area }, numel(C.nodes), 1);
    endif

    assert(numel(C.area) == numel(C.nodes), 'Invalid size of parameter <area>');
    t_area = cellfun(@(v) numel(v), C.area);

    assert(all(arrayfun(@(i) any(t_area(i) == [0, t_nodes(i)]), 1:numel(t_area))), 'The size of some area vector does not match the number of nodes');

    % Fields are now set. Some additionnal checks are performed
    lCmdType = HT_Cmd_GetType(C);

    if strcmpi(C.type, 'temperature')
      if ~C.force
        % If this is a temperature type command, the area resistance of the face is
        % check to be 0, otherwise an error is raised. If the resistance is not null
        % it means the face is not in direct contact with the nodes.
        t = cellfun(@(face) isempty(face.r) || (isnumeric(face.r) && all(abs(face.r) == eps('double'))), C.nodes);
        assert(all(t), sprintf('Cmd_Init::In command <%s>: Can not specify temperature on face <%s>. Thermal resistance between surface and nodes must be 0. Use option <force=true> to do it anyway.', ...
                                                          C.name, ...
                                                          HT_StrJoin(HT_GetStructField(C.nodes(~t), 'name'), ';')));
      endif

      % Make sure no area is specified 'none' if multiple identical nodes are found in a face.
      % In that case, SolveModel will compute the average temperature based on node area.
      if any(cellfun(@(v) isempty(v), C.area)) && any(lCmdType.type == [ HT_BTYPE_T_CONSTANT_NON_UNIFORM, HT_BTYPE_T_VARIABLE_NON_UNIFORM, HT_BTYPE_T_INTERPOLATION_ARRAY])
        error(sprintf('Cmd_Init::In command <%s>: Node areas must be specified when using NON_UNIFORM temperature boundary condition', C.name));
      endif
    endif
  endif
endfunction

