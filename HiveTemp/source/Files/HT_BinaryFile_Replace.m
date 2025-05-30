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
%  Copyright (c) 2025 AltRD-Emmanuel Ruffio
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
% 09/09/2024
% Replace data at <offset> of size <size> in file <file> by data stored in
% <data> at <dataoffset> of size <datasize>
% Save the result into target
%
% File, data, target may be a file handle or a string
function R = HT_BinaryFile_Replace(file, offset, size, data, dataOffset, dataSize, target)
unwind_protect
  lOutFile = 0;
  lSrcFile = 0;
  lDataFile = 0;
  R = 0;

  nOps = max([numel(offset) numel(size) numel(dataOffset) numel(dataSize)]);
  if isscalar(offset),      offset = repmat(offset, nOps, 1); endif
  if isscalar(size),        size = repmat(size, nOps, 1); endif
  if isscalar(dataOffset),  dataOffset = repmat(dataOffset, nOps, 1); endif
  if isscalar(dataSize),    dataSize = repmat(dataSize, nOps, 1); endif

  if ischar(file)
    lSrcFile = fopen(file, 'rb');
  else
    lSrcFile = file;
  endif

  if ischar(data)
    lDataFile = fopen(data, 'rb');
  else
    lDataFile = file;
  endif

  if ischar(target)
    lOutFile = fopen(target, 'wb');
  else
    lOutFile = file;
  endif

  assert(lOutFile > 0, 'Invalid target file');
  assert(lSrcFile > 0, 'Invalid source file');
  assert(lDataFile > 0, 'Invalid data file');

  % Go in position
  lStatus = fseek(lSrcFile, 0, SEEK_SET);
  assert(lStatus == 0, 'source fseek fails');
  % Read the data from the source file up to the first operation
  lData = fread(lSrcFile, offset(1));
  fwrite(lOutFile, lData); % and write to target file

  for i=1:(nOps-1)

    % Retrieve the data from the datafile
    lStatus = fseek(lDataFile, dataOffset(i), SEEK_SET);
    assert(lStatus == 0, 'data fseek fails');
    lData = fread(lDataFile, dataSize(i));
    assert(numel(lData) == dataSize(i));
    fwrite(lOutFile, lData);

    % Skip the block size specified by the user
    lStatus = fseek(lSrcFile, size(i), SEEK_CUR);
    assert(lStatus == 0, 'source fseek fails');

    % Read to the next offset operation
    lReadCount = offset(i+1)-offset(i)-size(i);
    assert(lReadCount >= 0);
    lData = fread(lSrcFile, lReadCount);
    assert(numel(lData) == lReadCount);
    fwrite(lOutFile, lData);
  endfor

  % Retrieve the data from the datafile
  lStatus = fseek(lDataFile, dataOffset(nOps), SEEK_SET);
  assert(lStatus == 0, 'data fseek fails');
  lData = fread(lDataFile, dataSize(nOps));
  assert(numel(lData) == dataSize(nOps));
  fwrite(lOutFile, lData);

  % Skip the block size specified by the user
  lStatus = fseek(lSrcFile, size(nOps), SEEK_CUR);
  assert(lStatus == 0, 'source fseek fails');

  lData = fread(lSrcFile);
  fwrite(lOutFile, lData);

unwind_protect_cleanup
  if ischar(file) && (lSrcFile > 0)
    fclose(lSrcFile);
  endif

  if ischar(target) && (lOutFile > 0)
    fclose(lOutFile);
  endif

  if ischar(data) && (lDataFile > 0)
    fclose(lDataFile);
  endif
end_unwind_protect
