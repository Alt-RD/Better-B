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
function [D titles header] = HT_ReadCsvFile(_file, varargin)
  if nargin < 2
    error('Missing input arguments');
  endif

  header = {};

  lParameters = struct( 'cacheFile', '', ...
                        'maxHeaderLine', 100, ...
                        'comment', '#', ...
                        'delimiter', ',', ...
                        'options', '', ...
                        'columns', [], ...
                        'forceReload', false, ...
                        'strict', true, ...   % Abort if error occur
                        'fillColumn', true, ... % If strict is false, fill the missing column with <fillValue>
                        'fillValue', NA, ...
                        'bufferLine', 1000, ...
                        'timeColumn', [], ... Specify the time column (necessary of subsampling if the format is not epoch)
                        'timeFormat', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                        'subSampling', 1);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("============ Read csv file ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'cacheFile',    '', @(v) ischar(v));
  lParameters = HT_CheckField(lParameters, 'options',    '', @(v) any(strcmpi({'', 'scan', 'reject'}, v)));
  lBackupParameters = lParameters;

  if ~isempty(lParameters.cacheFile) && ~lParameters.forceReload && isfile(lParameters.cacheFile)
    lFileData = load(lParameters.cacheFile);

    if isfield(lFileData, 'parameters') && isequaln(lFileData.parameters, lParameters)
##    if isfield(lFileData, 'parameters') && isfield(lFileData.parameters, 'columns') && isequal(lFileData.parameters, lParameters)
      D = lFileData.data;
      titles = lFileData.titles;
      header = lFileData.header;
      return;
    endif
  endif

  if ischar(lParameters.options) && isempty(lParameters.options)
    lParameters.options = {};
  elseif ischar(lParameters.options)
    lParameters.options = { lParameters.options };
  endif

  lRejectColumnMode = any(strcmpi(lParameters.options, 'reject'));
  lScanMode = any(strcmpi(lParameters.options, 'scan')) || ~isargout(1);

  unwind_protect
    if ischar(_file)
      lFileHandle = fopen(_file, 'r');
    else
      lFileHandle = _file;
      lParameters.maxHeaderLine = lParameters.bufferLine; % If a file handle is provided, the option maxHeaderLine is discarded
    endif

    lData = {};
    lDataTitles = {};

    lCurrentLine = 1;
    [lCurrentTitleFormat, lCurrentTitleLine, header, lInfos] = HT_CsvGuessFormat(lFileHandle, 'maxLineSearch', lParameters.maxHeaderLine, 'comment', lParameters.comment, 'delimiter', lParameters.delimiter);
    assert(~isempty(lCurrentTitleLine), sprintf('Could not find the title line in file <%s>', disp(_file)));
    lCurrentLine += lInfos.lineCount-1;

    while  ~feof(lFileHandle)
      % Find next title line
      [lTitleFormat, lTitleLine, ~, lInfos] = HT_CsvGuessFormat(lFileHandle, 'maxLineSearch', lParameters.bufferLine, 'comment', lParameters.comment, 'delimiter', lParameters.delimiter);

      lBlockLength = lInfos.skipRange(2) - lInfos.skipRange(1);
      if lBlockLength <= 0
        lCurrentTitleFormat = lTitleFormat;
        lCurrentTitleLine = lTitleLine;
        lCurrentLine += lInfos.lineCount-1;
        if feof(lFileHandle)
          break;
        endif

        fseek(lFileHandle, lInfos.blockRange(2), SEEK_SET());
        continue; % Means there were no data in the block
      endif

      % Before processing the title line, the block of data is processed
      fseek(lFileHandle, lInfos.skipRange(1), SEEK_SET());
      lDataBlock = fread(lFileHandle, lInfos.skipRange(2) - lInfos.skipRange(1), 'char=>char')';
      if ~lScanMode
        D = textscan(lDataBlock, strjoin(lCurrentTitleFormat), "CommentStyle", lParameters.comment, "Delimiter", lParameters.delimiter, "BufSize", numel(lDataBlock));
      else
        D = cell(1, numel(lCurrentTitleFormat));
      endif

      if ~isempty(lParameters.columns) % Select only required columns (if set)
        lColSelect = cellfun(@(v) find(strcmpi(lCurrentTitleLine, v)), lParameters.columns, 'UniformOutput', false);
        lColSelect(cellfun(@(v) isempty(v), lColSelect)) = [];
        lColSelect = cell2mat(lColSelect);

        if lRejectColumnMode
          D(lColSelect) = [];
          DTitles = lCurrentTitleLine;
          DTitles(lColSelect) = [];
        else
          D = D(lColSelect);
          DTitles = lCurrentTitleLine(lColSelect);
        endif
      else
        DTitles = lCurrentTitleLine;
      endif

      % Merge the new data with previous one
      if isempty(lData)
        lData = D;
        lDataTitles = DTitles;
      else
##        if isempty(lTitleLine) % Only data were found, not title line. It means the format is the same
##          lTitleLine = lCurrentTitleLine;
##          lTitleFormat = lCurrentTitleFormat;
##        endif

        [lData, lDataTitles, msg, newCol, ~] = HT_CsvMergeColumn(lData, lDataTitles, D, DTitles, struct('position', 'after'), lOptions);
        if ~isempty(msg), cellfun(@(v) disp(sprintf('Msg at line %d: %s', lCurrentLine, v)), msg); endif;
      endif

      % Skip the title line that was just read
      fseek(lFileHandle, lInfos.blockRange(2), SEEK_SET());

      if ~isempty(lTitleLine)
        lCurrentTitleFormat = lTitleFormat;
        lCurrentTitleLine = lTitleLine;
      endif
      lCurrentLine += lInfos.lineCount-1;
    endwhile

    % Apply subsampling if necessary
    if (lParameters.subSampling > 1)
      lIsTimeString = false;

      if ischar(lParameters.timeColumn)
        lParameters.timeColumn = find(strcmpi(lDataTitles, lParameters.timeColumn));
        assert(numel(lParameters.timeColumn) == 1, sprintf('Invalid time column name <%s>. Could not be found.', lParameters.timeColumn));
      endif

      if ~isempty(lParameters.timeColumn)
        lIsTimeString = ischar(D{lParameters.timeColumn}{1});
        if lIsTimeString
          t = datevec(D{lParameters.timeColumn}, lParameters.timeFormat);
          t = datenum(t);
          D(lParameters.timeColumn) = t;
        endif
      endif

      for i=1:numel(D)
        D{i} = HT_SubSampling(D{i}, lParameters.subSampling);
      endfor

      if lIsTimeString
        D(lParameters.timeColumn) = datestr(D{lParameters.timeColumn}, lParameters.timeFormat);
      endif
    endif

    if ~isempty(lParameters.cacheFile)
      [lCacheDir, lCacheName, lCacheExt] = fileparts(lParameters.cacheFile);

      lCacheFullDir = make_absolute_filename(lCacheDir);
      if ~isempty(lCacheFullDir) && ~isfolder(lCacheFullDir)
        lStatus = mkdir(lCacheFullDir);
        assert(islogical(lStatus) && lStatus, sprintf('Could not create cache directory <%s>', lCacheFullDir));
      endif

      parameters = lBackupParameters;
      data = lData;
      titles = lDataTitles;
      save('-binary', lParameters.cacheFile, 'data', 'titles', 'parameters', 'header');
    endif

    if lOptions.verbose, disp("Read done..."); endif;

    titles = lDataTitles;
    D = lData;

  unwind_protect_cleanup
    if ischar(_file) && (lFileHandle != 0)
      fclose(lFileHandle);
    endif
  end_unwind_protect

endfunction

function [titleLine blockRange state comments] = Int_FindTitleLine(_fileHandle, _parameters, _options, _state)
  titleLine = '';
  blockRange = [ftell(_fileHandle) NA];
  state = _state;
  comments = {};

  % Search for the column title line
  for i=1:_parameters.bufferLine
    tmp = fgetl(_fileHandle);
    state.currentLine++;

    if startswith(tmp, _parameters.comment)
      if isargout(4), comments = [comments; tmp]; endif
      continue;
    endif

    tmp = strtrim(tmp);
    if ~startsWith(tmp, _parameters.comment), break; endif;
  endfor

  if i != _parameters.bufferLine
    Int_Message(_options, sprintf('New column title line found at line %d', _state.currentLine));
  endif
endfunction

function Int_Message(options, str)
  if options.verbose
    disp(str);
  endif
endfunction

