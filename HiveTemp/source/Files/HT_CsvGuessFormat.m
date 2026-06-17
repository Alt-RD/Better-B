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
% F contains the format identifier (cell array): %f%q ...
% T contains the column titles
function [F T comments infos] = HT_CsvGuessFormat(_file, varargin)
  if nargin < 1
    error('Missing input arguments');
  endif

  F = [];
  T = [];
  comments = {};
  infos = struct('lineCount', 0, 'blockRange', [NA; NA], 'skipRange', [NA; NA]);

  lParameters = struct( 'maxLineSearch', 100, ...
                        'comment', '#', ...
                        'delimiter', ',', ...
                        'default', [] ...
                        );

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("== Guessing csv file format : "); disp('Checking input arguments'); endif;

  unwind_protect
    if ischar(_file)
      lFileHandle = fopen(_file, 'r');
    else
      lFileHandle = _file;
    endif

    lColumnTitle = [];
    lDataLine = [];
    lDataLineCount = 0;
    infos.blockRange(1) = ftell(lFileHandle);
    infos.skipRange(1) = infos.blockRange(1);

    % Find the title line (which contains the name of each column of data)
    for i=1:lParameters.maxLineSearch
      tmpFilePos = ftell(lFileHandle);
		  tmp = fgetl(lFileHandle);
      if feof(lFileHandle)
        break;
      endif

      infos.lineCount++;

      if startsWith(tmp, lParameters.comment)
        if isargout(3), comments = [comments; tmp]; endif;
        continue;
      endif

      tmp = strtrim(tmp);

      if startsWith(tmp, strcat(lParameters.comment, 'octave'))
        ind = index(tmp, '=');
        if (ind > 0) && (ind < numel(tmp))
          tmp = HT_StrTrimLeft(tmp(imd:end), "# ");
          if lOptions.verbose
            disp(sprintf("Field <octave> found: %s [...]", tmp(1:min(numel(tmp), 16))));
          endif
          [~, ~, ~, m, ~, ~, ~] = regexpi(tmp, '(%[a-z])');
          F = m;

          if isempty(F)
            disp('Invalid field');
          endif
        endif
      endif

      % First line of data ? it is column titles
      lIsTitleLine = Int_IsTitleLine(tmp, lParameters);
      %~startsWith(tmp, lParameters.comment) && ~isempty(tmp);
      if lIsTitleLine
        if lDataLineCount > 0 % Some data were found before the title line
          % Return the block of data before the current line
          infos.blockRange(2) = tmpFilePos;
          infos.skipRange(2) = tmpFilePos;
          break;
        elseif isempty(lColumnTitle)
          infos.skipRange(2) = tmpFilePos;
          lColumnTitle = tmp;
        else % Two successive title line were found
          infos.blockRange(2) = tmpFilePos;
          lDataLine = [];
          fseek(lFileHandle, tmpFilePos, SEEK_SET());
          break;
        endif
      else
        lDataLineCount++;
        if ~isempty(lColumnTitle)
          infos.blockRange(2) = tmpFilePos;
          lDataLine = tmp;
          fseek(lFileHandle, tmpFilePos, SEEK_SET());
          break;
        endif
      endif

      %infos.lineCount++;
    endfor

    if isna(infos.blockRange(2)), infos.blockRange(2) = ftell(lFileHandle); endif
    if isna(infos.skipRange(2)), infos.skipRange(2) = ftell(lFileHandle); endif

    if ~isempty(lColumnTitle)
      T = strsplit(HT_StrTrimLeft(lColumnTitle, "# "), lParameters.delimiter);
      T = cellfun(@(v) strtrim(v), T, 'uniformoutput', false);
      lColumnCount = numel(T);

      % No field <octave> was specified
      if isempty(F) && ~isempty(lDataLine)
        [~, ~, te, ~, ~, ~, ~] = regexp(lDataLine, '"(.*?)"');
        if !isempty(te)
          for k=1:numel(te)
            lDataLine(te{k}(1):te{k}(2)) = '0';
          endfor
        endif

        if !isempty(lParameters.delimiter)
          lDataFields = strsplit(lDataLine, lParameters.delimiter, 'collapsedelimiters', false);
        else
          lDataFields = strsplit(lDataLine, "[\\s\\t]+", 'collapsedelimiters', true, "delimitertype", "regularexpression")
        endif

        F = cellfun(@(v) Int_GetFormat(v), lDataFields, 'uniformoutput', false);
      endif

      if ~isempty(lDataLine) && (numel(F) != numel(T))
        error(sprintf('Inconsistent number of columns. Title line has %d columns. Data has %d.', numel(F), numel(T)));
      endif
    endif

    if lOptions.verbose
      if isempty(T) && (lDataLineCount > 0)
        disp(sprintf("Csv file: no title line found in the next %d lines, but %d data line found", lParameters.maxLineSearch, lDataLineCount));
      elseif isempty(T)
        disp(sprintf("Csv file: no title line found and no data found"));
      elseif isempty(F)
        disp(sprintf("Csv file contains %d columns but format could not be identified: No data line", lColumnCount));
      else
        disp(sprintf("Csv file contains %d columns with format <%s>", lColumnCount, strjoin(F)));
      endif
    endif

  unwind_protect_cleanup
    if ischar(_file) && (lFileHandle != 0)
      fclose(lFileHandle);
    endif
  end_unwind_protect

endfunction

function F = Int_GetFormat(_str)
  _str = strtrim(_str);
  if startsWith(_str, '"')
    F = '%q';
  elseif strcmpi(_str, "nan")
    F = '%f';
  else
    lValue = str2double(_str);
    if isnan(lValue)
      F = '%q';
    else
      F = '%f';
    endif
  endif
endfunction

function R = Int_IsTitleLine(_str, _parameters)
  % Retrieve the first element
  ind = strfind(_str, _parameters.delimiter);
  if isempty(ind)
    ind = numel(_str);
  else
    ind = ind(1);
  endif

  if isempty(_str)
    R = false;
    return;
  endif

  _str = _str(1:ind);
  if ((_str(1) >= '0') && (_str(1) <= '9')) || (_str(1) == '"')
    R = false;
  else
    R = true;
  endif
endfunction

