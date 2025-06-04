%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University
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
%
% Returns the temperature of a node or list of nodes from the matrix of temperature
% Input arguments:
% 1) <nodes> = [string or cell array of string] containing the node name or the list
%              of node names whose temperature must be extracted.
% 2) <Tmatrix> = [matrix dim Nnodes x Ntimes]
% 3) <Tnodes> = [cell array of string dim Nnodes]
% 4) <options>
%     .index = [vector of integer] If defined, returns only temperature whose
%              time index are defined in <.index>
%     .time = [vector] If defined, returns only temperature whose times are defined
%             in <.time>. The difference between <.index> and <.time> is the use
%             of interpolation. If a component could not be found in <.timeVector>
%             it is interpolated.
%     .timeVector = [vector] Column vector containing the time steps used by the
%             simulation
%     .merge = true/[false] If true, existing cell array in node list are expanded
%             first and all data are put in the same matrix. THe returned type is
%             a matrix.
%             If false, the returned type is a cell array whose element are matrices
%             corresponding to the nodes list
%
% Output arguments:
% Temperature matrix T (dim nNodes x nt)
function [T] = HT_Result_GetNodeTemperature(nodes, Tmatrix, Tnodes, varargin)
  HT_ImportConstants();

  assert(nargin >= 3, 'Missing input arguments');
  assert(ischar(nodes) || iscell(nodes), 'Invalid nodes');
  assert(isnumeric(Tmatrix) && (size(Tmatrix,1) == numel(Tnodes)), 'Invalid temperature matrix or node list');

  % Convert node name to cell array of name (of size 1)
  if ischar(nodes), nodes = { nodes }; endif

  prop = varargin(1:2:end);
  value = varargin(2:2:end);
  assert(numel(prop) == numel(value), 'Invalid properties');

  if HT_VAR_STRICT_INPUT_PARAMETERS
    lFieldNames = {'index', 'time', 'timevector', 'merge'};
    lFieldTest = prop;
    tf = ismember(lFieldTest, lFieldNames);
    assert(all(tf), sprintf('Invalid parameter name: %s', strjoin(cellfun(@(v) sprintf('%s',v), lFieldTest(tf), 'UniformOutput', false))));
    clear lFieldNames lFieldTest tf;
  endif

  lParams = struct( 'index', 1:size(Tmatrix,2), ...
                    'time', [], ....
                    'timevector', [], ...
                    'merge', true);

  for k=1:numel(prop)
    lParams = setfield(lParams, prop{k}, value{k});
  endfor

  if isempty(lParams.time), lParams.time = lParams.timevector; endif;

  if lParams.merge
    % If one element of <nodes> is a cell array ?
    if any(cellfun(@(v) iscell(v), nodes))
      nodes = HT_CellUnwrap(nodes);
      assert(all(cellfun(@(v) ischar(v), nodes)), 'Invalid node names');
    endif

    if ~isempty(lParams.time) && ~isempty(lParams.timevector)
      T = NaN(numel(nodes), numel(lParams.time));

      for k=1:numel(nodes)
        lMatchIndex = find(strcmpi(Tnodes, nodes{k}));
        assert(~isempty(lMatchIndex), sprintf('Node name <%s> could not be found in the node list', nodes{k}));
        assert(numel(lMatchIndex)==1, 'Node name was found multiple times');
        T(k,:) = interp1(lParams.timevector, Tmatrix(lMatchIndex,:), lParams.time, 'linear');
      endfor
    elseif ~isempty(lParams.index)
      T = NaN(numel(nodes), numel(lParams.index));

      for k=1:numel(nodes)
        lMatchIndex = find(strcmpi(Tnodes, nodes{k}));
        assert(~isempty(lMatchIndex), sprintf('Node name <%s> could not be found in the node list', nodes{k}));
        assert(numel(lMatchIndex)==1, 'Node name was found multiple times');
        T(k,:) = Tmatrix(lMatchIndex,lParams.index);
      endfor
    % Interpolate temperature based on times specified by the user
    else
      T = [];
    endif
  else % If lParams.merge
    T = {};
    for i=1:numel(nodes)
      T = [T; HT_Result_GetNodeTemperature(nodes{i}, Tmatrix, Tnodes, ...
                  'index', lParams.index, ...
                  'time', lParams.time, ...
                  'timevector', lParams.timevector, ...
                  'merge', true)];
    endfor
  endif

endfunction

