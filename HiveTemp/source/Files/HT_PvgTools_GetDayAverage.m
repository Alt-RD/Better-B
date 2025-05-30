%  This file is part of project HiveTemp.
%
%  Copyright (c) 2023 AltRD-Emmanuel Ruffio
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
% Returns the daily average data.
% Properties:
% 'group'= defines the base set over which the average is done
%          Only 'month' is implemented.
% 'mergeYear'= with group=month, it specifies if year must be merged together
function [t M Mmin Mmax Mstd] = HT_PvgTools_GetDayAverage(D, infos, varargin)

  props = varargin(1:2:end);
  values = varargin(2:2:end);

  lValidFields = {'group', 'merge'};
  assert(all(cellfun(@(v) any(strcmp(v, lValidFields)), props)), 'Invalid fields');

  lParams = struct('group', 'month', ...
                   'merge', false);

  for i=1:numel(props)
    lParams = setfield(lParams, props{i}, values{i});
  endfor

  t = [];
  M = [];

  if strcmp(lParams.group, 'month')
    % Number of values per days
    nPerDay = find(infos.datetime(:,3) > infos.datetime(1,3), 1)-1;
    lDateTime = double(infos.datetime);

    % Moyenne sur chaque mois
    if ~lParams.merge
      lMonthIndexVec = (lDateTime(:,1:2) - [lDateTime(1,1) 0]) * [12; 1];
      lMonthIndexList = lMonthIndexVec(1):lMonthIndexVec(end);
    else
      lMonthIndexVec = lDateTime(:,2);
      lMonthIndexList = unique(lMonthIndexVec);
    endif
%    lMonthIndexVec(lMonthIndexVec > 12) = [];

    lMonthCount = numel(lMonthIndexList);
    lColCount = columns(D);

    M = zeros(nPerDay*lMonthCount, lColCount);
    Mmin = [];
    Mmax = [];
    Mstd = [];

    if isargout(3), Mmin = zeros(nPerDay*lMonthCount, lColCount); endif
    if isargout(4), Mmax = zeros(nPerDay*lMonthCount, lColCount); endif
    if isargout(5), Mstd = zeros(nPerDay*lMonthCount, lColCount); endif

    for i=1:lMonthCount
      lMonthSet = find(lMonthIndexVec == lMonthIndexList(i));
      % Dcurrent contains only the data of the current month
      Dcurrent = D(lMonthSet,:);
      % Data are now splitted into days
      for k=1:lColCount
        Mt = Dcurrent(:,k);
        Mt = reshape(Mt, nPerDay, idivide(rows(Mt), int32(nPerDay)));
        lSet = (1+(i-1)*nPerDay):i*nPerDay;
        M(lSet,k) = mean(Mt,2);

        if isargout(3), Mmin(lSet,k) = min(Mt,[], 2); endif
        if isargout(4), Mmax(lSet,k) = max(Mt,[], 2); endif
        if isargout(5), Mstd(lSet,k) = std(Mt, 0, 2); endif
      endfor
    endfor

    t = lDateTime(1,[4 5]) * [1; 1/60]/24 + (0:(rows(M)-1))/24;

  elseif strcmp(lParams.group, 'year')
    % Number of values per year
    lYearChangeVec = 1;

    while true
      ind = find(infos.datetime(:,1) > infos.datetime(lYearChangeVec(end),1), 1);
      if isempty(ind), break; endif;

      lYearChangeVec = [lYearChangeVec; ind];
    endwhile

    lYearChangeVec = [lYearChangeVec; rows(infos.datetime)];
    lYearChangeVec = [lYearChangeVec(1:(end-1)), lYearChangeVec(2:end)-1];

    % Compute the number of data point for each years
    lYearSerieLength = diff(lYearChangeVec, 1, 2)+1;
    % Take the minimum (remove leap years)
    lSerieLength = min(lYearSerieLength);

    % Compute the average
    M = zeros(lSerieLength, rows(lYearChangeVec));
    for i=1:rows(lYearChangeVec)
      M(:,i) = D(lYearChangeVec(i,1) + (0:(lSerieLength-1)));
    endfor

    if isargout(3), Mmin = min(M, [], 2); endif
    if isargout(4), Mmax = max(M, [], 2); endif
    if isargout(5), Mstd = std(M, 0, 2); endif

    M = mean(M, 2);
    t = (0:(lSerieLength-1))/lSerieLength*365;
  else
    error('Unknown parameter');
  endif


endfunction
