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

  header = ''; % Not implemented

  lParameters = struct( 'cacheFile', '', ...
                        'maxHeaderLine', 100, ...
                        'comment', '#', ...
                        'delimiter', ',', ...
                        'columns', [], ...
                        'forceReload', false, ...
                        'subSampling', 1);

  lOptions = struct('skipFieldCheck', false, ...
                    'verbose', true);

  [lKeys lValues lParameters lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin);

  if lOptions.verbose, disp("============ Read csv file ================"); disp('Checking input arguments'); endif;

  lParameters = HT_CheckField(lParameters, 'cacheFile',    '', @(v) ischar(v));

  if ~isempty(lParameters.cacheFile) && ~lParameters.forceReload && isfile(lParameters.cacheFile)
    lFileData = load(lParameters.cacheFile);

    if isfield(lFileData, 'parameters') && isfield(lFileData.parameters, 'columns') && isequal(lFileData.parameters, lParameters)
      D = lFileData.data;
      titles = lFileData.titles;
      header = lFileData.header;
      return;
    endif
  endif

  unwind_protect
    if ischar(_file)
      lFileHandle = fopen(_file, 'r');
    endif

    % Find the title line (which contains the name of each column of data)
    for i=1:lParameters.maxHeaderLine
      tmp = strtrim(fgetl(lFileHandle));
      if ~startsWith(tmp, lParameters.comment) && ~isempty(tmp), break; endif;
      if ~isempty(tmp), lLine = tmp; endif
    endfor
    assert(i != lParameters.maxHeaderLine, 'Could not find title line');

    titles = strsplit(HT_StrTrimLeft(lLine, "# "), lParameters.delimiter);
    lColumnCount = numel(titles);

    if lOptions.verbose, disp(sprintf("Csv file contains %d columns", lColumnCount)); endif;

    lDataLineFormat = strjoin({'%s', repmat('%f', 1, lColumnCount-1)}, ',');

    D = textscan(lFileHandle, lDataLineFormat, "CommentStyle", lParameters.comment, "Delimiter", lParameters.delimiter);

    % Select only the specified columns if set
    if ~isempty(lParameters.columns)
      lIndexVec = cellfun(@(str) find(strcmpi(str, titles)), lParameters.columns, 'UniformOutput', false);
      lEmptyIndexVec = cellfun(@(v) isempty(v), lIndexVec);
      assert(all(~lEmptyIndexVec), sprintf('Could not find %d columns <%s>', numel(lEmptyIndexVec), strjoin(lParameters.columns(lEmptyIndexVec))));

      lIndexVec = cell2mat(lIndexVec);
      titles = titles(lIndexVec);
      D = D(1,lIndexVec);
    endif

    % Apply subsampling if necessary
    if lParameters.subSampling > 1
      for i=1:numel(D)
        D{1,i} = HT_SubSampling(D{1,i}, lParameters.subSampling);
      endfor
    endif

    if ~isempty(lParameters.cacheFile)
      [lCacheDir, lCacheName, lCacheExt] = fileparts(lParameters.cacheFile);

      lCacheFullDir = make_absolute_filename(lCacheDir);
      if ~isfolder(lCacheFullDir)
        lStatus = mkdir(lCacheFullDir);
        assert(islogical(lStatus) && lStatus, sprintf('Could not create cache directory <%s>', lCacheFullDir));
      endif

      parameters = lParameters;
      data = D;
      save(lParameters.cacheFile, 'data', 'titles', 'parameters', 'header');
    endif

    if lOptions.verbose, disp("Read done..."); endif;

  unwind_protect_cleanup
    if ischar(_file) && (lFileHandle != 0)
      fclose(lFileHandle);
    endif
  end_unwind_protect

endfunction

