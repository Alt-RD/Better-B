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
function [S binInfos] = HT_CsvBinAverage(D, T, varargin)
  if nargin < 3
    error('Missing input arguments');
  endif

  header = {};

  lParameters = struct( 'bin', 'month', ...
                        'unit', 'day', ...      % Char or hours
                        'maxHoleTime', [], ...  % Minutes
                        'columns', [], ...
                        'reject', [], ...
                        'timeColumn', '', ...
                        'timeType', '1970', ... % Used if 'timeColumn' is an epoch value
                        'time', int32(48), ...
                        'interpolation', 'linear', ... 'average'
                        'policy', 'loose', ...
                        'exportRawData', false, ...
                        'removeEmpty', true, ...
                        'statistics', @(V) std(V, 0, "omitnan")); % Function called to to a statistic analysis of each series

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("============ Get csv columns ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'bin',    'month', @(v) ischar(v) && strcmpi(v, 'month'));
  lParameters = HT_CheckField(lParameters, 'unit',   'day', @(v) ischar(v) && strcmpi(v, 'day'));
  lBackupParameters = lParameters;

  % Check columns names
  if isempty(lParameters.columns)
    lParameters.columns = T;
  endif

  if ischar(lParameters.columns)
    lParameters.columns = { lParameters.columns };
  endif

  if ~isempty(lParameters.reject)
    lRejectColumnList = lParameters.reject;
    if isnumeric(lRejectColumnList), lRejectColumnList = num2cell(lRejectColumnList); endif
    if ischar(lRejectColumnList), lRejectColumnList = { lRejectColumnList }; endif


    for i=1:numel(lRejectColumnList)
      if isnumeric(lRejectColumnList{i})
        lRejectColumnList{i} = T{lRejectColumnList{i}};
      endif
    endfor

    inds = cellfun(@(v) find(strcmpi(v, lParameters.columns), 1), lRejectColumnList, 'UniformOutput', false);
    assert(~any(cellfun(@(v) isempty(v), inds)));
    inds = cell2mat(inds);
    lParameters.columns(inds) = [];
  endif

  lColumnIndices = cellfun(@(v) HT_CellStrFilter(T, v), lParameters.columns, 'UniformOutput', false);
  lColumnIndicesNotFound = find(cellfun(@(v) isempty(v), lColumnIndices));
  assert(all(cellfun(@(v) (numel(v) <= 1), lColumnIndices)), 'Duplicate columns found after filter');

  assert(isempty(lColumnIndicesNotFound), sprintf('Sensors not found: %s', strjoin(lParameters.columns(lColumnIndicesNotFound))));
  assert(any(strcmpi(lParameters.interpolation, {'linear', 'average'})));

  if strcmpi(lParameters.interpolation, 'linear')
    lInterpolationFunc = @(x,y,t) interp1(x,y,t);
  elseif strcmpi(lParameters.interpolation, 'average')
    lInterpolationFunc = @(x,y,t) HT_Resample(x,y,t, 'interpolation', 'linear', 'type', 'point', 'average', true);
  else
    error('Invalid parameter');
  endif

  lColumnIndices = cell2mat(lColumnIndices);
  % Make sure all column are scalar values
  assert(all(cellfun(@(v) isvector(v) || isrow(v), D(lColumnIndices))), 'Some columns are not scalar values');

  % Check time column
  lTimeColIndex = lParameters.timeColumn;
  if ischar(lTimeColIndex)
    lTimeColIndex = find(strcmpi(T, lTimeColIndex));
    assert(numel(lTimeColIndex) == 1, sprintf('Invalid time column <%s>, could not be found', lParameters.timeColumn));
  endif

  % Split time column according to <bin> value
  lTimeCol = D{lTimeColIndex};
  % array of struct('datevec', 'range', 'index', 'name')
  [lTimeSplit lTimeEpoch] = HT_SplitTimeVector(lTimeCol, 'position', lParameters.bin, ...
                                            'boundaryPolicy', 'overlap', ...
                                            'removeEmpty', lParameters.removeEmpty, ...
                                            'epochType', lParameters.timeType);

  % For each bin, invalid time are detected (too large timestep)
  if ~isempty(lParameters.maxHoleTime)
    for i=1:numel(lTimeSplit)
      t = lTimeEpoch(lTimeSplit(i).index); % Index are sorted by time value (indices themselves are not necessarily sorted)

      lTimeSplit(i).invalidTime = (t(2:end)-t(1:end-1)) > (lParameters.maxHoleTime/(24*60)); % t is day number (datenum) and maxHoleTime is minute
    endfor
  endif

  % Build time vector
  lUnitWidthDay = Int_GetUnitWidth(lParameters); % Return the bin width in hours

  if isinteger(lParameters.time)
    lUnitTime = ((1:double(lParameters.time))-0.5)'/double(lParameters.time) * lUnitWidthDay;
  else
    lUnitTime = lParameters.time(:)/24;
    assert((numel(lUnitTime) > 1) && all((lUnitStep > 0) & (lUnitTime <= lUnitWidthDay)));
  endif

  %
  lIsLoose = strcmpi(lParameters.policy, 'loose');
  S = struct('name', cell(1:numel(lTimeSplit), 1:numel(lColumnIndices)), ...
             'data', [], ...
             't', [], ...
             'statistics', [], ...
             'rawData', []);

  for i=1:numel(lTimeSplit)
    lObj = lTimeSplit(i);
    if isempty(lObj.index), continue; endif % Empty bin ?

    tData = lTimeEpoch(lObj.index) - lObj.range(1);
    tUnitShift = 0:lUnitWidthDay:(lObj.range(2)-lObj.range(1) - lUnitWidthDay);

    % Extract the measurements of the current bin
##    M = cellfun(@(v) v(lObj.index), D(lColumnIndices), 'UniformOutput', false);
##    M = cell2mat(M);

    for k=1:numel(lColumnIndices)
      disp(k);
      y = D{lColumnIndices(k)}(lObj.index);
      y(lObj.invalidTime) = NA;
##      M = interp1(tData, y, lUnitTime + tUnitShift(:));
      M = lInterpolationFunc(tData, y, lUnitTime + tUnitShift);


      assert(lIsLoose || ~any(isna(y)));
      y = mean(M,2, "omitnan");
      S(i,k).name = sprintf("%s/%s", T{lColumnIndices(k)}, lObj.name);
      S(i,k).data = y;
      S(i,k).t = lUnitTime;
      if ~isempty(lParameters.statistics)
        S(i,k).statistics = arrayfun(@(i) lParameters.statistics(M(i,:)), 1:rows(M), 'UniformOutput', false);
        % Convert to vector if all values are numeric and scalar
        if all(cellfun(@(v) isscalar(v) && isnumeric(v), S(i,k).statistics))
          S(i,k).statistics = cell2mat(S(i,k).statistics)';
        endif
      endif

      if lParameters.exportRawData
        S(i,k).rawData = M;
      endif
    endfor
  endfor

  binInfos = lTimeSplit;

  if lOptions.verbose, disp("Bin average done..."); endif;

endfunction

% Return the unit duration (in days)
function W = Int_GetUnitWidth(_parameters)
  if ischar(_parameters.unit)
    if strcmpi(_parameters.unit, 'day')
      W = 1;
    else
      error('Not implemented');
    endif
  else
    assert(isnumeric(_parameters.unit) && isscalar(_parameters.unit));
    W = _parameters.unit/24;
  endif
endfunction


