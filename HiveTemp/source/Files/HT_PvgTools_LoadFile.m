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
% Conversion des fichiers de données https://re.jrc.ec.europa.eu/pvg_tools/en/
% de l'irradiance solaire.
% The matfile generated contains the following fields:
% .filename = the name of the file (if a string was specified as filename)
% .datetime = (uint dim Nx5) contains the Year, month, day, hour, min
% .t = (double dim Nx1) contains the time (in days)
% .D = (double dim Nx7) contains the data (beam radiation, diffuse radiation
%                       reflected radiation, elevation angle, air temperature,
%                       wind speed, flag
% .titles = (cell string dim 1x7) a short title for each column of D
% .descriptions = (cell string dim 1x7) a description for each column of D
% .stat = structure used to detect file change
% .comments = file header lines and file footer lines
%
% Retourne les données D et les titres des colonnes T d'un fichier
function [t D infos] = HT_PvgTools_LoadFile(_file)
  lFileNameSpecified = strcmp(class(_file), "char");
  lFile = [];

  unwind_protect

    % Try to load data from .mat file to avoid reloading all text file
    if lFileNameSpecified
      [lFileDir, lBaseFileName, lBaseFileExt] = fileparts(_file);
      lMatFilePath = strcat(lFileDir, '/', lBaseFileName, '.mat');

      % Compare the modification time
      [matstatus, materr] = stat(lMatFilePath);
      [lStat, err] = stat(_file);

      lReloadData = (materr != 0);
      % If ".mat" file does not exist, it is created
      if !materr % Check data are uptodate
        fileData = load("-mat", lMatFilePath);

        % Make sure all the following fields exist
        lFieldList = { 'descriptions', 'filename', 'D', 't', 'datetime', 'titles', 'stat', 'comments' };
        lFieldMissing = cellfun(@(v) ~any(strcmp(v, lFieldList)), fieldnames(fileData));
        if any(lFieldMissing)
          disp(sprintf('Some fields could not be found in file structure:', strjoin( lFieldList(lFieldMissing), ',' )));
          lReloadData = true;
        endif
        clear lFieldList lFieldMissing;

        lReloadData = lReloadData || (lStat.mtime != fileData.stat.mtime);
        lReloadData = lReloadData || (lStat.mtime != fileData.stat.mtime);
        if lReloadData, disp(sprintf("File <%s%s> has changed", lBaseFileName, lBaseFileExt)); endif;
      endif

      if ~lReloadData
        disp(sprintf('Reloading data from .mat file: <%s>', strcat(lBaseFileName, '.mat')));
        D = fileData.D;
        t = fileData.t;

        infos = struct( 'filename', fileData.filename, ...
                        'datetime', fileData.datetime, ...
                        'titles', { fileData.titles }, ...
                        'descriptions', { fileData.descriptions }, ...
                        'stat', fileData.stat, ...
                        'comments', { fileData.comments });

        return;
      endif
    endif

    % If data must be read from txt file
    if lFileNameSpecified % A filename
      lFile = fopen(_file, "rb");
      if lFile < 0, error(sprintf("Invalid file <%s>", _file)); endif;
    else
      lFile = _file;
    endif

    lFileOriginPos = ftell(lFile);
    lLinePos = lFileOriginPos;
    lComment = {};
    lLineCount = 1;

    % First step: looks for the first data line
    while ~feof(lFile)
      lLinePos = ftell(lFile);
      lLineStr = fgetl(lFile);

      % Try to extract date and time of the data line:
      % 20180101:0010,0.0,0.0,0.0,0.0,11.33,6.55,0.0
      [lValues lCount lErrMsg lPos] = sscanf(lLineStr, "%d:%d", 2);
      if (lCount == 2) && (lPos == 14)
        break;
      endif

      if lLineCount < 32
        lComment = [lComment; lLineStr];
      endif
      lLineCount++;
    endwhile

    if feof(lFile), error(sprintf('Could not find the first data line')); endif

    lDataPos = lLinePos;
    clear lLinePos lLineStr;

    if (fseek(lFile, lDataPos, SEEK_SET) != 0)
      error('Could not read the file');
    endif

    [C, ~, lErrMsg] = textscan(lFile, "%4u %2u %2u %2u %2u %f %f %f %f %f %f %f", ...
                            "Delimiter", ":,", ...
                            "ReturnOnError", true, ...
                            "CollectOutput", true);

    if ~isempty(lErrMsg), error(sprintf('Could not read the data : %s', lErrMsg)); endif

    if (C{1}(end,1) == 0), C{1}(end,:) = []; C{2}(end,:) = []; endif;

    lDataEnd = ftell(lFile);


    % Build the fields of the files
    % Time vector
    datetime = C{1};
    D = C{2};
    dt = datenum(double(datetime(2,:))) - datenum(double(datetime(1,:)));
    %t = arrayfun(@(i) datenum(double(datetime(i,:))), 1:rows(D));
    %t = t';
    t = dt * (1:rows(D))';
    filename = '';

    if lFileNameSpecified, filename = _file; endif;

    titles = {'Gb', ...
              'Gd', ...
              'Gr', ...
              'H_sun', ...
              'T2m', ...
              'WS10m', ...
              'Int' ...
              };

    descriptions = { 'Beam (direct) irradiance on the inclined plane (plane of the array) (W/m2)', ...
                    'Diffuse irradiance on the inclined plane (plane of the array) (W/m2)', ...
                    'Reflected irradiance on the inclined plane (plane of the array) (W/m2)', ...
                    'Sun height (degree)', ...
                    '2-m air temperature (degree Celsius)', ...
                    '10-m total wind speed (m/s)', ...
                    '1 means solar radiation values are reconstructed' ...
                    };

    comments = lComment;

    if (fseek(lFile, lDataEnd, SEEK_SET) != 0)
      error('Could not read the file');
    endif

    lLineCount = 1;

    % First step: looks for the first data line
    while ~feof(lFile) && (lLineCount < 32)
      lLineStr = fgetl(lFile);
      if ~isempty(lLineStr), comments = [comments; lLineStr]; endif;
      lLineCount++;
    endwhile

    if lFileNameSpecified
      % lFieldList = { 'descriptions', 'filename', 'D', 't', 'datetime', 'titles', 'stat', 'comments' };

      filename = _file;
      stat = stat(_file);
      save("-binary", "-mat", lMatFilePath, 'filename', 'datetime', 't', 'D', ...
                                            'titles', 'descriptions', 'stat', 'comments');
    else
      filename = '';
      stat = struct();
    endif

    infos = struct( 'filename', filename, ...
                    'datetime', datetime, ...
                    'titles', { titles }, ...
                    'descriptions', { descriptions }, ...
                    'stat', stat, ...
                    'comments', { comments });

  unwind_protect_cleanup
    if lFileNameSpecified && ~isempty(lFile) && (lFile >= 0), fclose(lFile); endif;
  end_unwind_protect
endfunction
