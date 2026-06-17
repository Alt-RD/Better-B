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
function [D T msg newCol missingCol] = HT_CsvMergeColumn(D, T, Dnew, Tnew, varargin)
  msg = {};

  lParameters = struct('position', 'after', ...
                       'titleFilter', @(s) s, ...
                       'emptyValue', NaN);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("== Merging csv column =="); endif;

  T = cellfun(@(v) lParameters.titleFilter(v), T, 'UniformOutput', false);
  Tnew = cellfun(@(v) lParameters.titleFilter(v), Tnew, 'UniformOutput', false);

  assert(numel(Tnew) == numel(unique(Tnew)), 'Duplicate column titles found in new data, once filtered');
  assert(isscalar(lParameters.emptyValue));

  lBuildEmptyArray = @(line, col) repmat(lParameters.emptyValue, line, col);

  % ==============================================================
  lColIndices = cellfun(@(v) find(strcmpi(T, v)), Tnew, 'UniformOutput', false);

  if !all(cellfun(@(v) numel(v) <= 1, lColIndices))
    msg = [msg; sprintf('Invalid title line <%s>', strjoin(Tnew, ';'))];
  endif

  lColNotFound = cellfun(@(v) isempty(v), lColIndices); % Index of new columns that do not exist in current data
  lColFound = ~lColNotFound;
  T = [T, Tnew(lColNotFound)];

##  lColIndices = cell2mat(lColIndices);

  lColNoData = 1:numel(D);% Index of columns for which there is no data in the current data
  lColNoData(cell2mat(lColIndices(lColFound))) = [];

  lColIndices(lColNotFound) = num2cell(numel(D)+1:numel(T));
  lColIndices = cell2mat(lColIndices);

  if ~isempty(D)
    lDataRowCount = rows(D{1});
  else
    lDataRowCount = 0;
  endif
  D(numel(D)+1:numel(T)) = lBuildEmptyArray(lDataRowCount,1); % Add empty columns to the data
##  D = [D, repmat({NA(lDataRowCount,1)}, 1, sum(lColNotFound))];

  % Add the new data into the old array
  D(lColIndices) = arrayfun(@(i) [D{lColIndices(i)}; Dnew{i}], 1:numel(lColIndices), 'UniformOutput', false);
##  for i=1:numel(lColIndices)
##    D{lColIndices(i)} = [D{lColIndices(i)}; D{i}];
##  endfor

  % Add empty data to columns that were not present in new data
  tmp = rows(Dnew{1});
  D(lColNoData) = arrayfun(@(i) [D{lColNoData(i)}; lBuildEmptyArray(tmp,1)], 1:numel(lColNoData), 'UniformOutput', false);
##  NA(tmp,1);
##  for i=1:numel(lColNoData)
##    D{lColNoData(i)} = [D{lColNoData(i)}; NA(tmp,1) ];
##  endfor

  newCol = lColNotFound;
  missingCol = lColNoData;

endfunction
