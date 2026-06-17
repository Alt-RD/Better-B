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
% Return a vector <u> where each element <i> is the integration of y(x) over
% tList(i) and tList(i+1) if 'type' is 'boundary'.
% <tList> may be a vector, a matrix, or a cellarray of matrix
% If option 'average' is set to true, <u> vector is divided by <w> to get the effective
%     (with NA/NAN removed) average.
% Output
% <w> is the effective width between tList(i) and tList(i+1). w(i) = tList(i+1)-tList(i)
%     If NAN or NA values are present, w will be smaller since invalid spaces are removed.
%     <w> has no NAN/NA values even if some are present in y, but <u> will be NA/NAN if
%     option 'omitNaN' is false.
%
function [u w] = HT_Resample(x,y,tList,varargin)
  lParameters = struct('interpolation', 'linear', ...
                       'type', 'boundary', ... % 'point'/'boundary'
                       'average', false, ...
                       'removeEmpty', false, ...
                       'omitNaN', false, ...
                       'strictMode', false);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  % ==============================================================
  if columns(x) > 1, x = x'; endif
  if columns(y) > 1, y = y'; endif

  assert(iscolumn(x) && iscolumn(y) && (numel(x) == numel(y)));
  assert(issorted(x));
  assert(any(strcmpi(lParameters.type, {'boundary', 'point'})));

  lIsBoundaryMode = strcmpi(lParameters.type, 'boundary');
  lIsCell = iscell(tList);
  lIsMatrix = isnumeric(tList) && (columns(tList) > 0);

  if ~lIsCell
    tList = mat2cell(tList, rows(tList), ones(columns(tList), 1));
  endif

  assert(all(cellfun(@(v) numel(v) > 1, tList)));

  if ~lIsBoundaryMode % Convert 'point' mode to 'boundary'
    p = cellfun(@(v) 0.5*(v(1:end-1)+v(2:end)), tList, 'UniformOutput', false);
    w = cellfun(@(v) diff(v), tList, 'UniformOutput', false);
    tList = arrayfun(@(i) [tList{i}(1)-w{i}(1)/2 ; p{i} ; tList{i}(end)+w{i}(end)/2], 1:numel(tList), 'UniformOutput', false);
  endif

  u = cell(1, numel(tList));

  for i=1:numel(tList)
    t = tList{i};

    if lParameters.strictMode
      assert((min(t) >= x(1)) && (max(t) <= x(end)));
    endif

    % Return the unique elements in t that are not in x.
    c = setdiff(t,x);
    if ~isempty(c)
      [xs, ixs] = sort([x; c]);
      ys = [y; interp1(x,y,c, lParameters.interpolation, NA)];
      ys = ys(ixs);
    endif
    ysInvalid = isna(ys) | isnan(ys);
    lIntegrated = 0.5*(ys(2:end)+ys(1:end-1)).*diff(xs);
    lIntegratedInvalid = ysInvalid(2:end) | ysInvalid(1:end-1); %isna(lIntegrated) | isnan(lIntegrated);

##      if lParameters.omitNaN
      lIntegrated(lIntegratedInvalid) = 0;
##      endif

    lIntegrated = [0; cumsum(lIntegrated)]; %cumtrapz(xs, ys);
##      lIntegrated = cumtrapz(xs, ys);

    [~, it] = ismember(t,xs);

    u{i} = lIntegrated(it(2:end)) - lIntegrated(it(1:end-1));

    if lParameters.average
      if any(lIntegratedInvalid)
##        if lParameters.omitNaN && any(lIntegratedInvalid)
        diffxs = diff(xs);
        diffxs(lIntegratedInvalid) = 0;
        cumdiffxs = cumsum([0; diffxs]);
        lWidth = cumdiffxs(it(2:end)) - cumdiffxs(it(1:end-1));
        w = lWidth;

        if ~lParameters.omitNaN
          lWidth(arrayfun(@(i) any(ysInvalid(it(i):it(i+1))), 1:numel(it)-1)) = NA;
        endif
        u{i} = u{i} ./ lWidth;
      else
        w = diff(xs(it));
        u{i} = u{i} ./ w;
      endif

    else
      if isargout(2), w = diff(xs(it)); endif
    endif
  endfor

  if lIsMatrix
    u = cell2mat(u);
  endif

endfunction



