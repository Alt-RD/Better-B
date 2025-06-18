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
function HT_WriteCsvFile(D, varargin)
  assert(nargin >= 2, 'Missing parameters');

  % TimeFormat may be "datevec", "epoch"
  % TimeColumn may be a column index (1x1) or an 1D/Matrix array containing the date
  lParameters = struct( 'header', '', ... %or "second" "epoch"
                        'footer', '', ...
                        'format', [], ...
                        'subSampling', 1, ...
                        'comment', '# ', ...
                        'linecount', [], ...
                        'titles', [], ...
                        'file', [], ...
                        'compress', false, ...
                        'separator', '', ...
                        'verbose', true);

  lOptions = struct('skipFieldCheck', false);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  nCol = columns(D);
  lParameters = HT_CheckField(lParameters, 'titles', [], {@(v) isempty(lParameters.titles) || numel(lParameters.titles) == nCol});

  if isempty(lParameters.titles) && ~ischar(lParameters.titles)
    lParameters.titles = arrayfun(@(v) sprintf('Col_%d', v), 1:nCol, 'UniformOutput', false);
  endif

  fLength = @(n) max(16,n+4); % Return the column title size with respect to the string length <n>
unwind_protect
  if ischar(lParameters.file)
    lFile = fopen(lParameters.file, 'w');
    assert(lFile != 0, sprintf('Invalid file <%s>', lParameters.file));
  else
    assert(lParameters.file > 0, 'Invalid file handle');
    lFile = lParameters.file;
  endif

  if ischar(lParameters.header) && ~isempty(lParameters.header)
    lParameters.header = { lParameters.header };
  endif

  if ~isempty(lParameters.header)
    lHeaderCells = HT_CellUnwrap(lParameters.header);
    for i=1:numel(lHeaderCells)
      assert(ischar(lHeaderCells{i}), sprintf('Invalid header data <%s>', lHeaderCells{i}));
      lHeaderCells{i} = strsplit(lHeaderCells{i}, "\n");
    endfor
    lHeaderCells = HT_CellUnwrap(lHeaderCells);
    cellfun(@(v) fputs(lFile, cstrcat(lParameters.comment, v, "\r\n")), lHeaderCells);
  endif

  % Compute the column width
  lColumnWidthVec = NA(nCol, 1);
  lFormatVec = cell(1, nCol);
  lFormatVec(1:numel(lParameters.format)) = lParameters.format;

  for i=1:nCol
    if iscell(D) && iscolumn(D)
      if iscell(D{i})
        assert(ischar(D{i}{1}));
        lLength = max(numel(D{i}{1}), numel(lParameters.titles{i}))+4;
        lFormatVec{i} = sprintf('%%%ds', lLength);
        lColumnWidthVec(i) = lLength;
      else
        assert(iscolumn(D{i}) && isnumeric(D{i}));
        [lFormatVec{i} lColumnWidthVec(i)] = INT_GetNumericFormat(lFormatVec{i}, max(abs(D{i})), numel(lParameters.titles{i}));
      endif
    elseif iscell(D)
      if ischar(D{1, i})
        lLength = max(numel(D{1, i}), numel(lParameters.titles{i}))+4;
        lFormatVec{i} = sprintf('%%%ds', lLength);
        lColumnWidthVec(i) = lLength;
      else
        assert(isnumeric(D{1,i}) && isscalar(D{1,i}));
        [lFormatVec{i} lColumnWidthVec(i)] = INT_GetNumericFormat(lFormatVec{i}, D{1,i}, numel(lParameters.titles{i}));
      endif
    else
      [lFormatVec{i} lColumnWidthVec(i)] = INT_GetNumericFormat(lFormatVec{i}, max(abs(D(:,i))), numel(lParameters.titles{i}));
    endif
  endfor

  % Write column titles
  if ~isempty(lParameters.titles)
    fprintf(lFile, '#');

    if ~lParameters.compress
      lTitleFormatVec = arrayfun(@(v) sprintf('%%%ds', v), lColumnWidthVec, 'UniformOutput', false);
      lTitleLine = strjoin(arrayfun(@(i) sprintf(lTitleFormatVec{i}, lParameters.titles{i}), 1:nCol, 'uniformoutput', false), '');
    else
      lTitleLine = strjoin(lParameters.titles, lParameters.separator);
    endif
    fputs(lFile, lTitleLine);
    fprintf(lFile, '\n');
  endif

  % Write data
##  lModelLine = strjoin(lFormatVec);
##  lModelLine = cstrcat(' ', lModelLine, '\n');
  lModelLine = strjoin(arrayfun(@(v) sprintf('%%%ds', v), lColumnWidthVec, 'UniformOutput', false), lParameters.separator);
  lModelLine = strcat(lModelLine, "\n");

  if iscell(D) && iscolumn(D) % Data variable is a cellarray containing vector/cell array
    lRowCount = min(cellfun(@(v) numel(v), D));
    if lParameters.verbose, disp(sprintf('Writing %d lines of data', lRowCount)); endif;

    for i=1:lRowCount
      lLineData = cellfun(@(v) INT_GetCellValue(v(i)), D, 'UniformOutput', false);
      fprintf(lFile, lModelLine, lLineData);

      if lParameters.verbose && mod(i, round(lRowCount/10)) == 0
        disp(sprintf('Writing line %d/%d', i, lRowCount));
      endif
    endfor
  else % Matrix or 2D cell array
    lRowCount = rows(D);
    lColCount = columns(D);
    if lParameters.compress
      for i=1:lRowCount
        lDataStr = arrayfun(@(k) sprintf(lFormatVec{k}, INT_GetCellValue(D(i,k))), 1:lColCount, 'UniformOutput', false);
        fputs(lFile, strjoin(lDataStr, lParameters.separator));
        fputs(lFile, "\n");

        if lParameters.verbose && mod(i, round(lRowCount/10)) == 0
          disp(sprintf('Writing line %d/%d', i, lRowCount));
        endif
      endfor
    else
      for i=1:lRowCount
        lDataStr = arrayfun(@(k) sprintf(lFormatVec{k}, INT_GetCellValue(D(i,k))), 1:lColCount, 'UniformOutput', false);
        fprintf(lFile, lModelLine, lDataStr{:});

        if lParameters.verbose && mod(i, round(lRowCount/10)) == 0
          disp(sprintf('Writing line %d/%d', i, lRowCount));
        endif
      endfor
    endif
  endif

  if ~ischar(lParameters.footer) && isempty(lParameters.footer)
    lParameters.footer = { repmat('=', 1, 80), ...
      sprintf('%d lines of data exported. %s', lRowCount, datestr(now)), ...
      repmat('=', 1, 80) };
  endif

  if ~isempty(lParameters.footer)
    cellfun(@(v) fprintf(lFile, cstrcat(lParameters.comment, v, '\n')), lParameters.footer);
  endif

unwind_protect_cleanup
  if (lFile > 0) && ischar(lParameters.file)
    fclose(lFile);
  endif
end_unwind_protect
endfunction

function V = INT_GetCellValue(V)
  if iscell(V)
    V = V{1};
  endif
endfunction

function [F columWidth] = INT_GetNumericFormat(defaultFormat, value, titleLength)
  if ~isempty(defaultFormat)
    F = defaultFormat;
    columWidth = max([numel(sprintf(defaultFormat, value)), titleLength+1]);
    return;
  endif

  columWidth = max([ceil(log10(max(1, value))) + 8, titleLength+1]);
  F = '%.4E';
endfunction

