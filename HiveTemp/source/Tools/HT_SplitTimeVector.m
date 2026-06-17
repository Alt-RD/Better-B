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
% Split a time vector (date string/epoch 1970/or datenum 0000)
% Return:
% <tinfos>, array of struct('datevec', 'range', 'index', 'name')
function [tinfos tepoch] = HT_SplitTimeVector(t, varargin)
  lParameters = struct('position', 'month', ...
                       'epochType', '1970', ... (seconds) or '0000' (day)
                       'timeFormat', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                       'boundaryPolicy', 'overlap', ...
                       'removeEmpty', false);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  % ==============================================================
  lIsCellStr = iscellstr(t);
  lIsEpochVec = isnumeric(t) && (iscolumn(t) || isrow(t));

  assert(lIsCellStr || lIsEpochVec, 'Invalid time vector');

  tepoch = Int_ConvertToEpoch(t, lParameters, lOptions);
  tinfos = Int_BuildSplitRanges(tepoch, lParameters, lOptions);
endfunction

function t = Int_ConvertToEpoch(t, lParameters, lOptions)
  if ~iscellstr(t)
    if strcmpi(lParameters.epochType, '1970')
      t = double(t)/86400 + datenum(1970,1,1,0,0,0);  %1970-01-01T00:00:00Z);
    endif
    return;
  endif

  t = datevec(t, lParameters.timeFormat);
  t = datenum(t);
  epochType = '0000';

endfunction

% Assume t is epoch 0000
% Returns S a struct array ('datevec', 'index', 'name', 'range')
function S = Int_BuildSplitRanges(t, lParameters, lOptions)
  lIsOverlapMode = strcmpi(lParameters.boundaryPolicy, 'overlap');
##  if lIsOverlapMode
    [tsorted tsortedind] = sort(t);
##  endif

  if strcmpi(lParameters.position, 'month')
    lMinTime = min(t);
    lMaxTime = max(t);

    % Round the min/max time to the previous month or next month
    lMinVec = datevec(lMinTime);
    lMaxVec = datevec(lMaxTime);
    lMaxVec = HT_AddToTimeVec(lMaxVec, 1, 'month');
    lMonthCount = (lMaxVec(1)-lMinVec(1))*12 + lMaxVec(2)-lMinVec(2);

    lEpochSplit = NA(lMonthCount+1, 1);
    lDateVecSplit = NA(lMonthCount+1, 6);
    x = [lMinVec([1 2]) 1 0 0 0];
    lEpochSplit(1) = datenum(x);
    lDateVecSplit(1,:) = x;

    for i=2:numel(lEpochSplit)
      x(2)++;
      if x(2) > 12
        x(1)++;
        x(2)-=12;
      endif
      lEpochSplit(i) = datenum(x);
      lDateVecSplit(i,:) = x;
    endfor

    if lIsOverlapMode
      lIndices = arrayfun(@(i) {find(tsorted>lEpochSplit(i), 1, 'first'), find(tsorted<lEpochSplit(i+1), 1, 'last') }, 1:lMonthCount, 'UniformOutput', false);
      lValuePerBin = cellfun(@(v) max(0, 1+v{2}-v{1}), lIndices);
      lIndices = cellfun(@(v) tsortedind(max(v{1}-1, 1):min(v{2}+1, numel(tsorted))), lIndices, 'UniformOutput', false);
      lIndices(lValuePerBin == 0) = {[]};
    else
      lIndices = arrayfun(@(i) find(tsorted>=lEpochSplit(i), 1, 'first'):find(tsorted<=lEpochSplit(i+1), 1, 'last'), 1:lMonthCount, 'UniformOutput', false);
      lValuePerBin = cellfun(@(v) numel(v), lIndices);
      lIndices = cellfun(@(v) tsortedind(v), lIndices, 'UniformOutput', false);
      lIndices(lValuePerBin == 0) = {[]};
##
##      lIndices = arrayfun(@(i) find(tsorted>=lEpochSplit(i)) & (tsorted<=lEpochSplit(i+1)), 1:lMonthCount, 'UniformOutput', false);
    endif

    S = struct('range', arrayfun(@(i) [lEpochSplit(i) lEpochSplit(i+1)], 1:lMonthCount, 'UniformOutput', false), ...
               'index', lIndices, ...
               'name', num2cell(datestr(lEpochSplit(1:lMonthCount), 28), 2)', ...
               'datevec', num2cell(lDateVecSplit(1:lMonthCount,:), 2)');
  else
    error('Not implemented');
  endif

  if lParameters.removeEmpty
    S(arrayfun(@(v) isempty(v.index), S)) = [];
  endif
endfunction

