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
% _file1 may be a file handle, a file name
% _ops a structure array: struct('title', '', 'description', '', 'function', @(i, D1, D2), 'source1', {}, 'source2', {})
function [data titles header] = HT_CsvOperation(_files1, _files2, _ops, varargin)
  if nargin < 3
    error('Missing input arguments');
  endif

  data = {};
  titles = {};
  header = ''; % Not implemented

  lParameters = struct( 'maxHeaderLine', 100, ...
                        'comment', '#', ...
                        'delimiter', ',', ...
                        'timeColumn1', [], ...
                        'timeColumn2', [], ...
                        'timeFormat1', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                        'timeFormat2', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                        'timeOutputFormat', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                        'timeOutputTitle', '', ...
                        'prefix1', '', ...
                        'maxHoleTime', 1200/86400, ... [day]
                        'check', true, ...
                        'forceReload', false, ...
                        'removeEmptyLines', true, ...
                        'cacheFile', '');

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("============ Read csv file ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'cacheFile',    '', @(v) ischar(v));

  if ~isempty(lParameters.cacheFile) && ~lParameters.forceReload && isfile(lParameters.cacheFile)
    lFileData = load(lParameters.cacheFile);

    if isfield(lFileData, 'parameters') && isfield(lFileData.parameters, 'columns') && isequal(lFileData.parameters, lParameters)
      error('not implemented');
      return;
    endif
  endif

  if ischar(_files1), _files1 = { _files1 }; endif
  if ischar(_files2), _files2 = { _files2 }; endif

  for i=1:numel(_ops)
    if ischar(_ops(i).source1)
      if isempty(_ops(i).source1)
        _ops(i).source1 = {};
      else
        _ops(i).source1 = { _ops(i).source1 };
      endif
    endif
    if ischar(_ops(i).source2)
      if isempty(_ops(i).source2)
        _ops(i).source2 = {};
      else
        _ops(i).source2 = { _ops(i).source2 };
      endif
    endif
  endfor

  unwind_protect
    lFileHandle = 0;
    lFileHandleSource = [];

    lSource1 = struct('D', {}, 'T', {}, 'filename', '');
    lSource2 = struct('D', {}, 'T', {}, 'filename', '');

    lRejectColFlag1 = arrayfun(@(v) any(strcmpi(v.options1, 'reject')), _ops);
    lRejectColFlag2 = arrayfun(@(v) any(strcmpi(v.options2, 'reject')), _ops);

    if any(lRejectColFlag1)
      if lOptions.verbose, disp(sprintf('Some operations reject columns. Scanning %d source files 1 to retrieve column titles', numel(_files1))); endif

      T = Int_GetAllColumnTitles(_files1, lParameters);

      if lOptions.verbose, disp(sprintf(' -> %d columns retrieved', numel(T))); endif

      if lParameters.check
        lRejectColumnList = HT_CellUnwrap({_ops(lRejectColFlag1).source1});
        lColumnExistList = cellfun(@(v) any(strcmpi(T, v)), lRejectColumnList);
        assert(all(lColumnExistList), sprintf('Some columns are rejected but do not exist: <%s>', strjoin(lRejectColumnList(~lColumnExistList), ';')));
        clear lRejectColumnList lColumnExistList;
      endif

      lRejectColFlag1 = find(lRejectColFlag1);

      for i=1:numel(lRejectColFlag1)
        _ops(lRejectColFlag1(i)).source1 = T(cellfun(@(v) ~any(strcmpi(_ops(lRejectColFlag1(i)).source1, v)), T));
      endfor

      clear T;
    endif

    if any(lRejectColFlag2)
      if lOptions.verbose, disp(sprintf('Some operations reject columns. Scanning %d source files 2 to retrieve column titles', numel(_files2))); endif

      T = Int_GetAllColumnTitles(_files2, lParameters);

      if lOptions.verbose, disp(sprintf(' -> %d columns retrieved', numel(T))); endif

      if lParameters.check
        lRejectColumnList = HT_CellUnwrap({_ops(lRejectColFlag2).source2});
        lColumnExistList = cellfun(@(v) any(strcmpi(T, v)), lRejectColumnList);
        assert(all(lColumnExistList), sprintf('Some columns are rejected but do not exist: <%s>', strjoin(lRejectColumnList(~lColumnExistList), ';')));
        clear lRejectColumnList lColumnExistList;
      endif

      lRejectColFlag2 = find(lRejectColFlag2);

      for i=1:numel(lRejectColFlag2)
        _ops(lRejectColFlag2(i)).source1 = T(cellfun(@(v) ~any(strcmpi(_ops(lRejectColFlag2(i)).source1, v)), T));
      endfor

      clear T;
    endif

    % Merge all column names
    lColumns1 = unique(HT_CellUnwrap({_ops.source1}));
    lColumns2 = unique(HT_CellUnwrap({_ops.source2}));

    % There is a special rule, if columns are specified in operation.source2 with no file in the source 2 list
    % it is assumed these columns are in source 1 list. It allows some special operations on columns.
    if isempty(_files2)
      lColumns1 = unique([lColumns1; lColumns2]);
      lColumns2 = [];
    endif

    if lParameters.check
      if numel(unique(lower(lColumns1))) != numel(lColumns1), warning(sprintf('Some columns differ only by lower/upper case')); endif
      if numel(unique(lower(lColumns2))) != numel(lColumns2), warning(sprintf('Some columns differ only by lower/upper case')); endif
    endif

##    lRejectColFlag1 = arrayfun(@(v) any(strcmpi(v.options1, 'reject')), _ops);
##    lRejectColFlag2 = arrayfun(@(v) any(strcmpi(v.options2, 'reject')), _ops);
##
##    % Merge all column names from source 1
##    lColumns1 = unique(HT_CellUnwrap({_ops(~lRejectColFlag1).source1}));
##    lColumns2 = unique(HT_CellUnwrap({_ops(~lRejectColFlag2).source2}));
##
##    lRemoveCol1 = unique(HT_CellUnwrap({_ops(lRejectColFlag1).source1}));
##    lRemoveCol2 = unique(HT_CellUnwrap({_ops(lRejectColFlag2).source2}));
##
##    % Column that must be removed but are requested by other operations are just removed
##    lRequiredCol1 = cellfun(@(v) any(strcmpi([lParameters.timeColumn1; lColumns1], v)), lRemoveCol1);
##    lRequiredCol2 = cellfun(@(v) any(strcmpi([lParameters.timeColumn2; lColumns2], v)), lRemoveCol2);
##    lRemoveCol1(lRequiredCol1) = [];
##    lRemoveCol2(lRequiredCol2) = [];

    % Read data from all files

    for i=1:numel(_files1)
      lFileHandleSource = _files1{i};

      if ischar(lFileHandleSource)
        lFileHandle = fopen(_files1{i}, 'r');
        assert(lFileHandle > 0, sprintf('Could not open the file <%s>',  lFileHandleSource));
      else
        lFileHandle = lFileHandleSource;
      endif

      [D T] = HT_ReadCsvFile(lFileHandle, struct( 'maxHeaderLine', lParameters.maxHeaderLine, ...
                                          'comment', lParameters.comment, ...
                                          'delimiter', lParameters.delimiter, ...
                                          'columns', {unique([lParameters.timeColumn1; lColumns1(:)])}, ...
                                          'timeColumn', lParameters.timeColumn1, ...
                                          'forceReload', lParameters.forceReload, ...
                                          'strict', false, ...
                                          'fillColumn', true, ... % Fill missing columns with NA
                                          'subSampling', 0));

      if ~isempty(lParameters.prefix1)
        if iscell(lParameters.prefix1)
          lPrefix = lParameters.prefix1{i};
        else
          lPrefix = lParameters.prefix1;
        endif

        lPrefix = strrep(lPrefix, '%(index)', num2str(i));
        if ischar(lFileHandleSource)
          lPrefix = strrep(lPrefix, '%(file)', num2str(i));
        else
          lPrefix = strrep(lPrefix, '%(file)', sprintf('handle%d',i));
        endif

        lIsDataColumn = ~strcmpi(lParameters.timeColumn1, T);
        T(lIsDataColumn) = cellfun(@(v) strcat(lPrefix, '.', v), T(lIsDataColumn), 'UniformOutput', false);
      endif

      fclose(lFileHandle);
      lFileHandle = 0;
      lFileHandleSource = [];

      lFileName = '<handle>';
      if ischar(_files1{i}), lFileName = _files1{i}; endif
      lSource1 = [lSource1; struct('D', {D}, 'T', {T}, 'filename', lFileName)];
    endfor

    for i=1:numel(_files2)
      lFileHandleSource = _files2{i};

      if ischar(lFileHandleSource)
        lFileHandle = fopen(_files2{i}, 'r');
        assert(lFileHandle > 0, sprintf('Could not open the file <%s>', lFileHandleSource));
      else
        lFileHandle = lFileHandleSource;
      endif

      [D T] = HT_ReadCsvFile(lFileHandle, struct( 'maxHeaderLine', lParameters.maxHeaderLine, ...
                                          'comment', lParameters.comment, ...
                                          'delimiter', lParameters.delimiter, ...
                                          'columns', {unique([lParameters.timeColumn2; lColumns2(:)])}, ...
                                          'timeColumn', lParameters.timeColumn2, ...
                                          'forceReload', lParameters.forceReload, ...
                                          'strict', false, ...
                                          'fillColumn', true, ... % Fill missing columns with NA
                                          'subSampling', 0));

      fclose(lFileHandle);
      lFileHandle = 0;
      lFileHandleSource = [];

      lFileName = '<handle>';
      if ischar(_files2{i}), lFileName = _files2{i}; endif
      lSource2 = [lSource2; struct('D', {D}, 'T', {T}, 'filename', lFileName)];
    endfor

    % Merge all data of source 1
    if ~isempty(lSource1)
      data1 = {};
      titles1 = {};

      for i=1:numel(lSource1)
        [data1, titles1] = HT_CsvMergeColumn(data1, titles1, lSource1(i).D, lSource1(i).T);
      endfor

      lTimeCol1 = lParameters.timeColumn1;
      if ischar(lParameters.timeColumn1)
        lTimeCol1 = find(strcmpi(titles1, lParameters.timeColumn1));
        assert(numel(lTimeCol1) == 1, sprintf('Invalid time column 2 <%s>, could not be found', lParameters.timeColumn1));
      endif

      time1 = data1{lTimeCol1};
      if iscell(time1)
        time1 = datevec(time1, lParameters.timeFormat1);
      endif
      if columns(time1) == 6
        time1 = datenum(time1);
      endif
      assert(columns(time1) == 1, sprintf('Invalid time column 1'));
      if ~issorted(time1, 'ascend')
        lDiff = time1(2:end)-time1(1:end-1);
        lDiff = find(lDiff < 0);
        lDiff = lDiff(1:min(4, numel(lDiff)));
        lInvalidTimeList = data1{lTimeCol1}(lDiff);
        lIterFunc = {@arrayfun, @cellfun}{1+iscell(lInvalidTimeList)};
        error(sprintf('Time1 is not sorted, around %d time samples: Example <%s>', numel(lDiff), strjoin(lInvalidTimeList, ';')));
        clear lDiff lInvalidTimeList lIterFunc;
      endif
    endif

    % Merge all data of source 1
    data2 = {};
    titles2 = {};
    time2 = {};

    if ~isempty(lSource2)
      data2 = {};
      titles2 = {};

      for i=1:numel(lSource2)
        [data2, titles2] = HT_CsvMergeColumn(data2, titles2, lSource2(i).D, lSource2(i).T);
      endfor

      lTimeCol2 = lParameters.timeColumn2;
      if ischar(lParameters.timeColumn2)
        lTimeCol2 = find(strcmpi(titles2, lParameters.timeColumn2));
        assert(numel(lTimeCol2) == 1, sprintf('Invalid time column 2 <%s>, could not be found', lParameters.timeColumn2));
      endif

      time2 = data2{lTimeCol2};
      data2(lTimeCol2) = [];
      titles2(lTimeCol2) = [];

      if iscell(time2)
        time2 = datevec(time2, lParameters.timeFormat1);
      endif
      if columns(time2) == 6
        time2 = datenum(time2);
      endif
      assert(columns(time2) == 1, sprintf('Invalid time column 2'));
      assert(issorted(time2, 'ascend'), 'Invalid time vector 2. Not sorted');

      % Synchronise source 2 with source 1
      [lInter1 lInterInd1] = HT_CsvSplitWithHoleTime(time1, lParameters.maxHoleTime);
      [lInter2 lInterInd2] = HT_CsvSplitWithHoleTime(time2, lParameters.maxHoleTime);

      assert(~isempty(lInter1) || ~isempty(lInter2), "No valid time range were found. Check maxHoleTimeValue");

      [lInterInd1, lInterInd2] = Int_IntersectSets(lInter1, lInterInd1, lInter2, lInterInd2);

      lSubSetTime1 = NA(numel(time1), 1);
      lSubSetTime1(lInterInd1) = time1(lInterInd1);

      data2 = cellfun(@(v) interp1(time2, v, lSubSetTime1), data2, 'UniformOutput', false);
    elseif isempty(_files2)
      % If no source 2 is specified, it is implicitly set to source 1
      data2 = data1;
      titles2 = titles1;
      time2 = time1;
    endif

##lOperations = struct( 'title', 'TestDiff', ...
##                      'description', '', ...
##                      'function', @(i, D1, D2) D1 - D2, ...
##                      'source1', {'34AluminumCover.Wall.T(T190@0.48)', 'AirExt.T(T112@3.4f)'}, ...
##                      'source2', '');

##    T = HT_CellUnwrap({_ops.outputTitle});

    % Apply operations
    for i=1:numel(_ops)
      lOp = _ops(i);

      if lOptions.verbose, disp(sprintf("Csv operation: processing operation <%s>", lOp.name)); endif;

      lData1 = [];
      lData2 = [];

      lInd1 = cellfun(@(v) find(strcmpi(titles1, v), 1), lOp.source1, 'UniformOutput', false);
      lInd2 = cellfun(@(v) find(strcmpi(titles2, v), 1), lOp.source2, 'UniformOutput', false);

      lColNotFound1 = cellfun(@(v) isempty(v), lInd1);
      lColNotFound2 = cellfun(@(v) isempty(v), lInd2);

      assert(all(~lColNotFound1), sprintf('Could not find source columns <%s> in source 1', strjoin(_ops(i).source1(lColNotFound1))));
      assert(all(~lColNotFound2), sprintf('Could not find source columns <%s> in source 2', strjoin(_ops(i).source2(lColNotFound2))));

      lInd1 = [lInd1{:}];
      lInd2 = [lInd2{:}];

      lData1 = [data1{lInd1}];
      lData2 = [data2{lInd2}];

      R = _ops(i).function(i, lData1, lData2);
      lOutputColumnNames = lOp.outputTitle;

      if ischar(lOp.outputTitle)
        lOutputColumnNames = {lOp.outputTitle};
      elseif is_function_handle(lOutputColumnNames) % Specify a name model
        lParams = struct('sourceName1', {lOp.source1}, 'sourceName2', {lOp.source2}, 'operationName', lOp.name);
        lOutputColumnNames = arrayfun(@(k) lOutputColumnNames(k, lParams), 1:columns(R), 'UniformOutput', false);
      elseif iscell(lOutputColumnNames)
        lOutputColumnNames = lOp.outputTitle{k};
      endif

      assert(columns(R) == numel(lOutputColumnNames), sprintf('Output matrix has %d columns. It does not match the number of outputs <%d>', columns(R), numel(lOp.outputTitle)));

      data = [data, num2cell(R, 1)];
      titles = [titles, lOutputColumnNames];
    endfor

    % Add the time
    if ~strcmpi(lParameters.timeOutputFormat, 'epoch')
      time1 = arrayfun(@(v) strftime('%Y/%m/%dT%H:%M:%S.000', gmtime(v)), time1, 'UniformOutput', false);
##      time1 = cellstr(datestr(time1, lParameters.timeOutputFormat));
    endif

    if ~isempty(lParameters.timeOutputTitle)
      time1Title = lParameters.timeOutputTitle;
    else
      time1Title = lParameters.timeColumn1;
    endif

    if lParameters.removeEmptyLines && ~isempty(data)
      lRemoveFlag = true(numel(data{1}), 1);
      for i=1:numel(data)
        lRemoveFlag = lRemoveFlag & (isna(data{i}) | isnan(data{i}));
      endfor

      time1(lRemoveFlag,:) = [];
      data = cellfun(@(v) v(~lRemoveFlag), data, 'UniformOutput', false);
    endif

    data = [{time1}, data];
    titles = [time1Title, titles];

    if lOptions.verbose, disp("Csv operation done..."); endif;

  unwind_protect_cleanup
    if ischar(lFileHandleSource) && (lFileHandle > 0)
      fclose(lFileHandle);
      lFileHandle = 0;
      lFileHandleSource = [];
    endif
  end_unwind_protect

endfunction

% Return the indices of vector elements <tset> that exist in at least of set of <trangeset>
function [setsInd1 setsInd2] = Int_IntersectSets(sets1, setsInd1, sets2, setsInd2)
  for i=1:numel(sets1)
    ind = Int_Intersect(sets1{i}, sets2);
    setsInd1{i} = setsInd1{i}(ind);
  endfor

  for i=1:numel(sets2)
    ind = Int_Intersect(sets2{i}, sets1);
    setsInd2{i} = setsInd2{i}(ind);
  endfor

  setsInd1(cellfun(@(v) isempty(v) || (numel(v) == 1), setsInd1)) = [];
  setsInd2(cellfun(@(v) isempty(v) || (numel(v) == 1), setsInd2)) = [];

  setsInd1 = cell2mat(setsInd1');
  setsInd2 = cell2mat(setsInd2');
endfunction


% Return the indices of vector elements <tset> that exist in at least of set of <trangeset>
function ind = Int_Intersect(tset, trangeSet)
  lValid = false(numel(tset), 1);
  for i=1:numel(trangeSet)
    lValid = lValid | ((tset >= trangeSet{i}(1)) & (tset <= trangeSet{i}(end)));
  endfor

  ind = lValid;
endfunction

% Return the total list of column titles found
function T = Int_GetAllColumnTitles(_fileList, _parameters)
  T = {};
  lFileHandle = 0;

  unwind_protect

  for i=1:numel(_fileList)
    if ischar(_fileList{i})
      lFileHandle = fopen(_fileList{i}, 'r');
      assert(lFileHandle > 0, sprintf('Could not open the file <%s>', _fileList{i}));
    else
      lFileHandle = _fileList{i};
    endif

    [~, tmp] = HT_ReadCsvFile(lFileHandle, struct( 'maxHeaderLine', _parameters.maxHeaderLine, ...
                                        'comment', _parameters.comment, ...
                                        'delimiter', _parameters.delimiter, ...
                                        'columns', {{}}, ... Read all columns
                                        'forceReload', _parameters.forceReload, ...
                                        'strict', false, ...
                                        'subSampling', 0));
    if ischar(_fileList{i})
      fclose(lFileHandle);
      lFileHandle = 0;
    endif

    T = [T, tmp];
  endfor

  T = unique(T);

  unwind_protect_cleanup
    if lFileHandle > 0
      fclose(lFileHandle);
      lFileHandle = 0;
    endif
  end_unwind_protect
endfunction

% Build the output column name
% _infos fields are: 'columnIndex', 'source1Name', 'source2Name', 'operationName'
function STR = Int_BuildOutputColumnName(_titleStr, _infos)
  [~,~,~, ~, t] = regexp(_titleStr, '\$\(([\.\w]+)\)');
  t = HT_CellUnwrap(t);
##  lFieldList = cellfun(@(v) strcat('$(', v, ')'), fieldnames(_infos), 'UniformOutput', false);
  lFieldList = fieldnames(_infos);
  lFieldCheckList = cellfun(@(v) any(strcmp(lFieldList, v)), t);
  assert(all(lFieldCheckList), sprintf('Invalid output column fields <%s>', strjoin(t(lFieldCheckList), '; ')));

  for i=1:numel(t)
    _titleStr = strrep(_titleStr, strcat('$(', t{i}, ')'), getfield(_infos, t{i}));
  endfor

  STR = _titleStr;
endfunction
