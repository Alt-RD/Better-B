%  This file is part of project HiveTemp.
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
% Returns temperatures related to a face
% Input arguments:
% 1) <face object> = [object or array of type Face]
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
%     .operation = ["average"], "array"
%     .merge = [logical] [false] if true and if face is an array, all
%               data are merged into a single value (average of average for
%               "average" operation) or array.
%
%
% Output arguments:
% Temperature matrix T (dim nNodes x nt)
function [T] = HT_Result_GetFaceTemperature(face, Tmatrix, Tnodes, varargin)
  HT_ImportConstants();

  assert(nargin >= 3, 'Missing input arguments');
  assert(HT_CheckType(face, 'face') || ...
          (iscell(face) && all(cellfun(RT_CheckType(v, 'face'), face))) || ...
          (isstruct(face) && all(arrayfun(RT_CheckType(v, 'face'), face))), 'Invalid face object or face array');
  assert(isnumeric(Tmatrix) && (rows(Tmatrix) == numel(Tnodes)), 'Invalid temperature matrix or node list');

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

  % Retrieve the temperature of each nodes for each faces
  T = cell(numel(face), 1);
  S = NA(numel(face), 1); % Store the total area of each faces

  for i=1:numel(T)
    lFaceObject = Int_GetObject(face(i));

    Tvec = HT_Result_GetNodeTemperature(lFaceObject.nodes, ...
                                        Tmatrix, ...
                                        Tnodes, ...
                                        'index', lExtractNodeParams.index, ...
                                        'time', lExtractNodeParams.time, ...
                                        'timevector', lExtractNodeParams.timevector);

    if strcmpi(lParams.operation, 'average')
      Snodes = HT_Face_GetAbsoluteData(lFaceObject, 'nodesArea');
      S(i) = sum(Snodes);
      T{i} = (Snodes' / S(i)) * Tvec;
    elseif strcmpi(lParams.operation, 'array')
      T{i} = Tvec;
    else
      error(sprintf('Operation <%s> not implemented', lParams.operation));
    endif
  endfor

  if lParams.merge
    if strcmpi(lParams.operation, 'array')
      T = cell2mat(T(:));
    elseif strcmpi(lParams.operation, 'average')
      if numel(T) > 1
        T = cell2mat(T(:));
        T = (T .* S) / sum(S);
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

