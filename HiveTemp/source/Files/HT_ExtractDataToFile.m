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
% 2nd version of buildBinjrFile.
% This version allows to extract only a subpart of sensor
function [R titles] = HT_ExtractDataToFile(_file, varargin)
  if nargin < 1
    error('Missing input arguments');
  endif

  % TimeFormat may be "datevec", "epoch"
  % TimeColumn may be a column index (1x1) or an 1D/Matrix array containing the date
  lParameters = struct( 'timeFormat', 'timestamp', ... %or "second" "epoch"
                        'timeColumn', [], ...
                        'timeFlag', [], ...
                        'data', [], ...
                        'sensorNames', [], ...   % Sensor names corresponding to <data>
                        'extractList', [], ...  % [struct] .name .title .digit .remove .calibration: List of sensor names to extract
                        'subSampling', 1,
                        'output', 'cell', ... or 'matrix'
                        'header', [], ...
                        'footer', [], ...
                        'separator', '', ...
                        'compress', false, ...
                        'verbose', true);

  lOptions = struct('skipFieldCheck', false);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lParameters.verbose, disp("============ Extracting data ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'timeFormat',    'timestamp', {});
  lParameters = HT_CheckField(lParameters, 'timecolumn',    [], {@(v) isempty(v) || iscolumn(v) || columns(v) == 6});
  lParameters = HT_CheckField(lParameters, 'data',          [],    {'exist'});
  lParameters = HT_CheckField(lParameters, 'sensorNames',   [],   {'exist', @(v) columns(lParameters.data) == numel(v)});
  lParameters = HT_CheckField(lParameters, 'extractList',   [],   {'exist', @(v) isstruct(v)});
  lParameters = HT_CheckField(lParameters, 'subSampling',   1,    {@(v) isnumeric(v)});
  lParameters = HT_CheckField(lParameters, 'output',        [],   {@(v) any(strcmpi(v, {'cell', 'matrix'}))});

  lParameters.sensorNames = HT_StrTrim(lParameters.sensorNames, " '\"");

  assert(!strcmpi(lParameters.timeFormat, 'timestamp') || strcmpi(lParameters.output, 'cell'), 'Only cell output is supported for time stamp column');

  if lParameters.verbose, disp('Building time vector...'); endif;

  if isscalar(lParameters.timeColumn) % If a scalar is provided, it is assumed it is the column index
    if iscell(lParameters.data)
      lTimeColumn = lParameters.data{lParameters.timeColumn};
    else
      lTimeColumn = lParameters.data(:,lParameters.timeColumn);
    endif
  else
    lTimeColumn = lParameters.timeColumn;
  endif

  % Convert the time column into the appropriate format
  if ~isempty(lTimeColumn)
    lTimeColumn = INT_ConvertTimeVector(lTimeColumn, lParameters);
    if ~isempty(lParameters.timeFlag)
      lTimeColumn = lTimeColumn(lParameters.timeFlag);
    endif
  endif

  unwind_protect
    lFileHandle = 0;

    if ischar(_file)
      lFileHandle = fopen(_file, 'w');
    elseif ~iscell(_file)
      lFileHandle = _file;
    endif

    nCol = 0;
    nRow = 0;

    if lParameters.verbose, disp('Exporting data...'); endif;

    % Find the index of each sensor in <extractList>
    lDataIndex = arrayfun(@(v) find(strcmpi(v.name, lParameters.sensorNames), 1), lParameters.extractList, 'uniformOutput', false);
    lDataIndexEmpty = cellfun(@(v) isempty(v), lDataIndex);
    for i=find(lDataIndexEmpty)
      disp(sprintf('Could not find sensor name <%s>', lParameters.extractList(i).name));
    endfor
    lDataIndex(lDataIndexEmpty) = [];
    lDataIndex = cell2mat(lDataIndex);
    lValidSensor = find(~lDataIndexEmpty);

    R = cell(1, numel(lValidSensor));

    for i=1:numel(lValidSensor)
      ind = lDataIndex(i);
      v = lParameters.extractList(lValidSensor(i));

      if iscell(lParameters.data{ind}) % Column of cells, it is rejected
        if lParameters.verbose, disp(sprintf('Data column %d <%s> is rejected: contains cell array ', ind, v.title)); endif;
        lDataIndex = NA;
        continue;
      elseif columns(lParameters.data{ind}) > 1 % Column does not contain scalar values ?
        if lParameters.verbose, disp(sprintf('Data column %d <%s> is rejected: not scalar values', i, v.title)); endif;
        lDataIndex = NA;
        continue;
      endif

      if lParameters.verbose
        disp(sprintf('Exporting data column %d <sensor=%s> to column %d <columntitle=%s>', ind, v.name, i, v.title));
      endif

      if iscell(lParameters.data)
        R{i} = lParameters.data{ind}(:);
      else
        R{i} = lParameters.data(:,ind);
      endif

      if ~isempty(lParameters.timeFlag)
        R{i} = R{i}(lParameters.timeFlag);
      endif

      if ~isempty(v.calibration)
        R{i} = v.calibration(R{i});
      endif

      if ~isempty(v.remove)
        lRemoveFlag = false(size(R{i}));
        for k=1:numel(v.remove)
          if is_function_handle(v.remove{k})
            lRemoveFlag = lRemoveFlag | v.remove{k}(R{i});
          elseif ischar(v.remove{k}) && strcmpi(v.remove{k}, 'spike')
##            lRateVariation = (R{i}(1:(end-2)) - 2*R{i}(2:(end-1)) + R{i}(3:end)) > std(R{i});
            lOutliers = abs(R{i} - movmedian(R{i},5)) > 0.5*std(R{i}(~isnan(R{i}) & ~isna(R{i}) ));
            lRemoveFlag = lRemoveFlag | lOutliers;
          elseif all(isfloat(v.remove{k})) & numel(lParameters.timeColumn) == numel(R{i})
            lRemoveFlag = lRemoveFlag | (lParameters.timeColumn >= v.remove{k}(1) & lParameters.timeColumn <= v.remove{k}(2));
          else
            assert(all(isinteger(v.remove{k})), 'Invalid field <remove>');
            lRemoveFlag(v.remove{k}(1): v.remove{k}(2)) = true;
          endif
        endfor

        R{i}(lRemoveFlag) = NA;
      endif

      R{i} = HT_SubSampling(R{i}, lParameters.subSampling);
    endfor

    % Remove invalid columns
    lRemoveColumn = isna(lDataIndex);
    R(lRemoveColumn) = [];
    lDataIndexEmpty(lRemoveColumn) = true;
    lValidSensor = find(~lDataIndexEmpty);
    lDataIndex(lRemoveColumn) = [];

    titles = arrayfun(@(v) v.title, lParameters.extractList(lValidSensor), 'UniformOutput', false);

    if ~isempty(lTimeColumn) % If the timestamp is provided in a separate vector one column must be added
      titles = [lParameters.timeFormat, titles];
      R = [cell(1,1), R];
      R{1,1} = lTimeColumn;
    endif

    if strcmpi(lParameters.output, 'cell')
      nCol = columns(R);
      nRow = median(cellfun(@(v) numel(v), R));

      M = cell(nRow, nCol);

      for i=1:columns(R)
        if ~iscell(R{i})
          M(:,i) = mat2cell(R{i}, ones(numel(R{i}), 1));
        else
          M(:,i) = R{i};
        endif
      endfor
    else % Matrix output
      nCol = columns(R);
      nRow = median(cellfun(@(v) numel(v), R));

      M = NA(nRow, nCol);
      for i=1:columns(R)
        M(:,i) = R{i};
      endfor
    endif

    R = M;

    if lFileHandle != 0
      lUserFormatVec = {lParameters.extractList.digit};

      lCharFormatVec = cellfun(@(v) ischar(v), lUserFormatVec);
      lDigitFormatVec = cellfun(@(v) isnumeric(v), lUserFormatVec);
      lColumnFormatVec = cell(1, numel(lUserFormatVec));
      lColumnFormatVec(lCharFormatVec) = lUserFormatVec(lCharFormatVec);
      lColumnFormatVec(lDigitFormatVec) = cellfun(@(v) sprintf('%%.%df', v), lUserFormatVec(lDigitFormatVec), 'UniformOutput', false);

      if ~isempty(lTimeColumn), lColumnFormatVec = ['%s', lColumnFormatVec]; endif;
      HT_WriteCsvFile(R,'file', lFileHandle, ...
                        'titles', titles, ...
                        'format', lColumnFormatVec, ...
                        'header', lParameters.header, ...
                        'compress', lParameters.compress, ...
                        'separator', ',',...
                        'footer', lParameters.footer);
    endif

    if strcmpi(lParameters.output, 'cell')
      R = [cell(1, columns(R)); R];
      R(1,:) = titles;
    endif

    if lParameters.verbose, disp("Exporting done..."); endif;

  unwind_protect_cleanup
    if ischar(_file)
      fclose(lFileHandle);
    endif
  end_unwind_protect

endfunction

function S = INT_setOptionStructure(_defaultStructure, _varargs)
  props = _varargs(1:2:end);
  values = _varargs(2:2:end);

  lOptionList = fieldnames(_defaultStructure);
  assert(all(cellfun(@(v) any(strcmpi(v, lOptionList)), props)), 'Invalid input argument');

  S = _defaultStructure;

  for i=1:numel(props)
    ind = find(strcmpi(props{i}, lOptionList));
    assert( ~isempty(ind), 'Invalid input field');
    S = setfield(S, lOptionList{ind}, values{i});
  endfor

endfunction

function lTimeColumn = INT_ConvertTimeVector(lTimeColumn, lParameters)
  if lParameters.verbose
    disp('Converting time vector to appropriate string format...');
  endif

  if (lParameters.subSampling != 1)
    assert(iscolumn(lTimeColumn), "Subsampling is not supported with timevector");
    lTimeColumn = HT_SubSampling(lTimeColumn, lParameters.subSampling);
  endif

  lValidFlag = ~any(isna(lTimeColumn) | isnan(lTimeColumn), 2);

  % Set the time vector in the correct format
  if strcmpi(lParameters.timeFormat, 'timestamp')
    if columns(lTimeColumn) == 6 % Assumed to be datetime vector
      lTimeColumn(lValidFlag) = cellstr(datestr(lTimeColumn(lValidFlag), 31));
      lTimeColumn(~lValidFlag) = NA;
    elseif isvector(lTimeColumn) % Assumed to be epoch values
      lColumnStr = cellstr(datestr(lTimeColumn(lValidFlag), 31));
      lTimeColumn = repmat({''}, numel(lTimeColumn), 1);
      lTimeColumn(lValidFlag) = lColumnStr;
    else
      error('Unknown time vector format');
    endif
  elseif strcmpi(lParameters.timeFormat, 'second')
    if columns(lTimeColumn) == 6 % Assumed to be datetime vector
      lTimeColumn(lValidFlag) = datenum(lTimeColumn(lValidFlag));
      lTimeColumn(~lValidFlag) = NA;
      lTimeColumn -= lTimeColumn(1);
      lTimeColumn *= 86400; % Convert from day to second
    elseif isvector(lTimeColumn) % Assumed to be epoch values
      lTimeColumn -= lTimeColumn(1);
      lTimeColumn *= 86400; % Convert from day to second
    else
      error('Unknown time vector format');
    endif
  elseif strcmpi(lParameters.timeFormat, 'epoch')
    if columns(lTimeColumn) == 6 % Assumed to be datetime vector
      lTimeColumn(lValidFlag) = datenum(lTimeColumn(lValidFlag));
      lTimeColumn(~lValidFlag) = NA;
    elseif isvector(lTimeColumn) % Assumed to be epoch values
    else
      error('Unknown time vector format');
    endif
  endif
endfunction

