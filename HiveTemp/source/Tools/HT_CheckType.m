%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University, AltRD-Emmanuel Ruffio
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
% This functions check the type of specified structure.
% Returns <true> if <obj> is compatible with type <typename>
%
function T = HT_CheckType(objList, typename)
  lTypeList = {};

  if ischar(typename)
    lTypeList = { typename };
  elseif iscell(typename)
    lTypeList = typename;
  else
    error('<typename> must be string or a cell array of string');
  endif

  T = false(numel(objList), 1);

  if isstruct(objList)
    for i=1:numel(objList)
      obj = objList(i);
      if isfield(obj, '__type__') && any(strcmpi(obj.__type__, lTypeList))
        T(i) = true;
      endif
    endfor
  elseif iscell(objList)
    for i=1:numel(objList)
      obj = objList{i};
      if isstruct(obj) && isfield(obj, '__type__') && any(strcmpi(obj.__type__, lTypeList))
        T(i) = true;
      endif
    endfor
  else
    T = false;
  endif

endfunction

