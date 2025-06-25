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
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2023-2025: CoActions-AltRD-Emmanuel Ruffio
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
function [lKeys lValues lParams lOptions] = HT_ReadVarargin(lDefaultParams, lDefaultOptions, args, type)
  assert(nargin >= 3, 'Missing parameters');
  if nargin < 4, type = 'struct'; endif
  assert(any(strcmpi(type, {'struct', 'struct array'})));

  if isempty(lDefaultParams), lDefaultParams = struct(); endif
  if isempty(lDefaultOptions), lDefaultOptions = struct(); endif

  lKeys = {};
  lValues = {};
  lParams = lDefaultParams;
  lOptions = lDefaultOptions;

  if isempty(args)
    return;
  endif

  if isstruct(args{1})        % If first parameter is a struct, parameters are loaded from it
    lParams = INT_MergeStruct(lParams, args{1});

    if (mod(numel(args), 2) == 0)
      assert(isstruct(args{end}), 'Invalid arguments <args>. End param must be a struct');

      lOptions = args{end};
      lKeys = args(2:2:(end-1));
      lValues = args(3:2:(end-1));
    else
      lKeys = args(2:2:end);
      lValues = args(3:2:end);
    endif
  else
    if (mod(numel(args), 2) == 1)
      assert(isstruct(args{end}), 'Invalid arguments <args>. End param must be a struct');

      lOptions = INT_MergeStruct(lOptions, args{end});
      lKeys = args(1:2:(end-1));
      lValues = args(2:2:(end-1));
    else
      lKeys = args(1:2:end);
      lValues = args(2:2:end);
    endif
  endif

  % If a default structure was specified, make sure all properties lKeys
  % are found in that structure
  if ~isempty(lParams)
    lFieldNames = fieldnames(lParams);
    lValidKeys = cellfun(@(v) any(strcmpi(v, lFieldNames)), lKeys);
    assert(all(lValidKeys), sprintf('Some keys are invalid: %s', strjoin(lKeys(~lValidKeys))));
  endif

  if isargout(3)
    lFieldNames = fieldnames(lParams);

    % Add to the <lParams> structure all parameters specified
    for i=1:numel(lKeys)
      lParams = setfield(lParams, lKeys{i}, lValues{i});
      if ~any(strcmp(lFieldNames, lKeys{i}))
        lFieldNames = [lFieldNames; lKeys{i}];
      endif
    endfor

    if strcmpi(type, 'struct array')
      lFieldNames = fieldnames(lParams);

      lCells = cell(2*numel(lFieldNames), 1);
      lCells(1:2:end) = lFieldNames;
      lCells(2:2:end) = cellfun(@(v) getfield(lParams, v), lFieldNames, 'UniformOutput', false);
      lParams = struct(lCells{:});
    endif
  endif
endfunction

function obj1 = INT_MergeStruct(obj1, obj2)
  if isempty(obj1)
    obj1 = obj2;
    return;
  endif

  lField2 = fieldnames(obj2);
  for i=1:numel(lField2);
    assert(isfield(obj1, lField2{i}));
    obj1 = setfield(obj1, lField2{i}, getfield(obj2, lField2{i}));
  endfor
endfunction
