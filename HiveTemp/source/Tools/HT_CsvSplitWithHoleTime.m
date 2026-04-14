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
% Split the vector of time based on time difference
% If the time difference between 2 samples exceeds <maxHoleTime>, the vector
% is split.
% SETCELL is a cellarray of time value. Cell content is sorted and cells are sorted
% INDCELL is a cellarray containing the index of each value of SELCELL in the original
%         time vector.
function [SETCELL INDCELL] = HT_CsvSplitWithHoleTime(tvec, _maxHoleTime)
  if columns(tvec) == 6
    % Convert to datenum
    tvec = datenum(tvec);
  endif
  assert(columns(tvec) == 1, sprintf('Invalid time vector size. Number of columns is %d instead of 1', columns(tvec)));

  isort = [];
  if !issorted(tvec)
    [tvec isort] = sort(tvec, 'ascend');
  endif

  lDifft = tvec(2:end) - tvec(1:end-1);
  lSplit = find(lDifft > _maxHoleTime);

  if ~isempty(lSplit)
    INDCELL = cell(numel(lSplit)+1, 1);
    INDCELL{1} = 1:lSplit(1);

    for i=2:numel(lSplit)
      INDCELL{i} = (lSplit(i-1)+1):lSplit(i);
    endfor
    INDCELL{end} = (lSplit(end)+1):numel(tvec);
    INDCELL(cellfun(@(v) numel(v) == 1, INDCELL)) = []; % Remove sets than contain only one element
  else
    INDCELL = { 1:numel(tvec) };
  endif

  if isargout(1)
    SETCELL = cellfun(@(v) tvec(v), INDCELL, 'UniformOutput', false);
  endif

  if isargout(2)
    if ~isempty(isort)
      INDCELL = cellfun(@(v) isort(v), INDCELL, 'UniformOutput', false);
    endif
  endif

endfunction
