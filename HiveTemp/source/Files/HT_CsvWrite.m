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
function S = HT_CsvWrite(_file, _data, _titles, varargin)
  [lFormatList lFormatListTitles] = Int_BuildFormatList(_data, _titles);

  _parameters = struct( 'overwrite', false, ...
                        'headerFields', struct('key', {}, 'value', {}), ...
                        'delimiter', ',', ...
                        'comment', '#' ...
                        );

  _options = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [~,~, _parameters, _options] = HT_ReadVarargin(_parameters, _options, varargin);

  if _options.verbose, disp("============ Writing csv file ================"); endif;

  S = [];

 unwind_protect
   lFileHandle = 0;

  if ischar(_file)
    lWriteMode = {'a', 'w'}{1+_parameters.overwrite};
    lFileHandle = fopen(_file, lWriteMode);
    assert(lFileHandle > 0, sprintf('Could not open the file <%s>', _file));
  else
    lFileHandle = _file;
  endif

  lColumnCount = numel(_titles);
  lColumnSize = cellfun(@(v) numel(v), _data);
  lColumnSize = unique(lColumnSize);
  assert(numel(lColumnSize) == 1, sprintf('Invalid column length. They should have the same length. Got <%s>', sprintf('%d;', lColumnSize)));

  fprintf(lFileHandle, "%c ============================================================\n", _parameters.comment);
  fprintf(lFileHandle, "%c File created by HT_CsvWrite function\n", _parameters.comment);
  fprintf(lFileHandle, "%c Date=%s\n", _parameters.comment, datestr(now()));
  fprintf(lFileHandle, "%c Column count=%d\n", _parameters.comment, numel(_titles));
  fprintf(lFileHandle, "%c Row count=%d\n", _parameters.comment, lColumnSize);

  for i=1:numel(_parameters.headerFields)
    fprintf(lFileHandle, "%c %s=%s\n", _parameters.comment, _parameters.headerFields(i).key, _parameters.headerFields(i).value);
  endfor
  fprintf(lFileHandle, "%c ============================================================\n", _parameters.comment);

  lLine = arrayfun(@(k) sprintf(lFormatListTitles{k}, _titles{k}), 1:numel(_titles), 'UniformOutput', false);
  fprintf(lFileHandle, strcat(strjoin(lLine, _parameters.delimiter), "\n"));

  % Convert cell array of string to matrix of string to speed up write operations
  lStrColumn = find(cellfun(@(v) endsWith(v, 's'), lFormatList));
  for i=1:numel(lStrColumn)
    _data{lStrColumn(i)} = strvcat(_data{lStrColumn(i)});
  endfor

  for i=1:lColumnSize
    lLine = arrayfun(@(k) sprintf(lFormatList{k}, _data{k}(i,:)), 1:numel(_titles), 'UniformOutput', false);
    fprintf(lFileHandle, strcat(strjoin(lLine, _parameters.delimiter), "\n"));

    if (mod(i,100) == 0) && _options.verbose
      disp(sprintf("Writing data %d/%d", i, lColumnSize));
    endif
  endfor

  fprintf(lFileHandle, "%c ============================================================\n", _parameters.comment);

  if ischar(_file)
    fclose(lFileHandle);
    lFileHandle = 0;
  endif

  unwind_protect_cleanup
    if ischar(_file) && (lFileHandle > 0)
      fclose(lFileHandle);
      lFileHandle = 0;
    endif
  end_unwind_protect

endfunction

function [F Ftitles] = Int_BuildFormatList(_data, _titles)
  F = cell(numel(_data), 1);
  W = zeros(numel(_data), 1);

  for i=1:numel(_data)
    if iscell(_data{i})
      assert(all(cellfun(@(v) ischar(v), _data{i})));
      W(i) = max(cellfun(@(v) numel(v), _data{i}));
      F{i} = 's';
    elseif isnumeric(_data{i})
      t = _data{i};
      t = t(isfinite(t));
      if isempty(t)
        t = 0;
      else
        t = ceil(max(log10(1+abs(t))));
      endif
      t += 4;   % 4 means '.' and 2 decimals and potential sign
      W(i) = t;
      F{i} = '.2f';
    else
      error(sprintf('Invalid type at column %d', i));
    endif
  endfor

  lTitleWidthList = cellfun(@(v) numel(v), _titles) + 4;
  lColumnWidth = max([lTitleWidthList(:), W], [], 2);

  F = arrayfun(@(k) sprintf('%%%d%s', lColumnWidth(k), F{k}), 1:numel(_data), 'UniformOutput', false);
  Ftitles = arrayfun(@(k) sprintf('%%%ds', lColumnWidth(k)), 1:numel(_data), 'UniformOutput', false);
endfunction
