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
%  Copyright (c) 2023 AltRD-Emmanuel Ruffio
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
% Returns temperature of a volume
% Input arguments:
% 1) <volume object> = [object or array of type Face]
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
%             simulation, that correspond to the matrix columns of <Tmatrix>
%     .operation = ["average"], "array"
%     .merge = [logical] [false] if true and if volume is an array, all
%               data are merged into a single value (average of average for
%               "average" operation) or array.
%
%
% Output arguments:
% Temperature matrix T (dim nNodes x nt)
function [T infos] = HT_Result_GetVolumeTemperature(volume, Tmatrix, Tnodes, varargin)
  HT_ImportConstants();

  assert(nargin >= 3, 'Missing input arguments');
  assert(HT_CheckType(volume, 'volume') || ...
          (iscell(volume) && all(cellfun(@(v) HT_CheckType(v, 'volume'), volume))) || ...
          (isstruct(volume) && all(arrayfun(@(v) HT_CheckType(v, 'volume'), volume))), 'Invalid volume object');
  assert(isnumeric(Tmatrix) && (rows(Tmatrix) == numel(Tnodes)), 'Invalid temperature matrix or node list');

  T = [];
  infos = struct('V', []);

  prop = varargin(1:2:end);
  value = varargin(2:2:end);
  assert(numel(prop) == numel(value), 'Invalid properties');

  if HT_VAR_STRICT_INPUT_PARAMETERS
    lFieldNames = {'index', 'time', 'timevector', 'merge', 'operation'};
    lFieldTest = prop;
    tf = ismember(lFieldTest, lFieldNames);
    assert(all(tf), sprintf('Invalid parameter name: %s', strjoin(cellfun(@(v) sprintf('%s',v), lFieldTest(tf), 'UniformOutput', false))));
    clear lFieldNames lFieldTest tf;
  endif

  lParams = struct();
  for i=1:numel(prop)
    lParams = setfield(lParams, prop{i}, value{i});
  endfor

  lParams = HT_CheckField(lParams, 'operation',       'average', @(v) any(strcmpi(v, {"average", "array"})));
  lParams = HT_CheckField(lParams, 'merge',           'false', @(v) islogical(v));

  % Copy fields from varargin to struct lParams
  lExtractNodeParams = struct('index', [], 'time', [], 'timevector', []);
  lCopyField = fieldnames(lParams);

  for i=1:numel(prop)
    if any(strcmpi(prop{i}, lCopyField)), lExtractNodeParams = setfield(lExtractNodeParams, prop{i}, value{i}); endif;
  endfor

  % Retrieve the temperature of each nodes for each volume
  T = cell(numel(volume), 1);
  V = NA(numel(volume), 1); % Store the total area of each volume

  for i=1:numel(T)
    if iscell(volume) && iscell(volume{i})
      [T lInfos] = HT_Result_GetNodeTemperature(volume{i}, ...
                                                Tmatrix, ...
                                                Tnodes, ...
                                                'index', lExtractNodeParams.index, ...
                                                'time', lExtractNodeParams.time, ...
                                                'timevector', lExtractNodeParams.timevector, ...
                                                'operation', 'average');
      V(i) = lInfos.V;
      T{i} = T;
    else
      lVolumeObject = Int_GetObject(volume(i));

      Tvec = HT_Result_GetNodeTemperature(lVolumeObject.nodes, ...
                                          Tmatrix, ...
                                          Tnodes, ...
                                          'index', lExtractNodeParams.index, ...
                                          'time', lExtractNodeParams.time, ...
                                          'timevector', lExtractNodeParams.timevector);

      V(i) = sum(lVolumeObject.nodesVolume);
      T{i} = sum((lVolumeObject.nodesVolume / V(i)) .* Tvec, 1);
    endif

  endfor

  if lParams.merge
    if strcmpi(lParams.operation, 'array')
      T = cell2mat(T(:));
      infos.V = V;
    elseif strcmpi(lParams.operation, 'average')
      infos.V = sum(V);

      if numel(T) > 1
        T = cell2mat(T(:));
        T = sum((T .* V) / infos.V, 1);
      else
        T = T{1};
      endif
    else
      error('<Merge> operation is not impemented in this configuration');
    endif
  endif

endfunction

function obj = Int_GetObject(obj)
  if iscell(obj)
    obj = obj{1};
  endif
endfunction

