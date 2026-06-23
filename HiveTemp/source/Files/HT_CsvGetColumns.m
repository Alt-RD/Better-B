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
function [S t] = HT_ReadCsvFile(D, T, varargin)
  if nargin < 3
    error('Missing input arguments');
  endif

  header = {};

  lParameters = struct( 'cacheFile', '', ...
                        'timeColumn', '', ...
                        'columns', '', ...
                        'forceReload', false, ...
                        'subSampling', 1);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("============ Get csv columns ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'cacheFile',    '', @(v) ischar(v));
  lParameters = HT_CheckField(lParameters, 'columns',      {}, @(v) ischar(v) || iscellstr(v));
  lBackupParameters = lParameters;

  if ~isempty(lParameters.cacheFile) && ~lParameters.forceReload && isfile(lParameters.cacheFile)
    lFileData = load(lParameters.cacheFile);

    if isfield(lFileData, 'parameters') && isequaln(lFileData.parameters, lParameters)
      S = lFileData.S;
      t = lFileData.t;
      return;
    endif
  endif

  if ischar(lParameters.columns)
    lParameters.columns = { lParameters.columns };
  endif

  lColumnIndices = cellfun(@(v) HT_CellStrFilter(T, v), lParameters.columns, 'UniformOutput', false);
  lColumnIndicesNotFound = find(cellfun(@(v) isempty(v), lColumnIndices));

  assert(isempty(lColumnIndicesNotFound), sprintf('Sensors not found: %s', strjoin(lParameters.columns(lColumnIndicesNotFound))));

  lColumnIndices = cell2mat(lColumnIndices);

  S = struct('title', T(lColumnIndices), 'data', D(lColumnIndices));

  % Load the time column vector
  ind = strcmp(T, lParameters.timeColumn);
  assert(~isempty(ind), 'Invalid time column title');
  t = D{ind};

  if lOptions.verbose, disp("Get columns done..."); endif;

endfunction



