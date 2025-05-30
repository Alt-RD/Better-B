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
function S = HT_StrJoin(list, delimiter, maxsize, varargin)
  if nargin < 3, maxsize = Inf; endif;

  props = varargin(1:2:end);
  values = varargin(2:2:end);

  lParams = struct('doubleFormat', '%g', ...
                   'integerFormat', '%d');

  for i=1:numel(props)
    if strcmpi(props{i}, 'doubleFormat')
      lParams.doubleFormat = values{i};
    elseif strcmpi(props{i}, 'integerFormat')
      lParams.integerFormat = values{i};
    endif
  endfor

  lPrefix = '';

  if numel(list) > maxsize
    lPrefix = sprintf('+...[%d elts]', numel(list)-maxsize);
  endif

  n = min(numel(list), maxsize);

  if ~iscell(list), list = arrayfun(@(v) v, list, 'uniformOutput', false); endif

  for i=1:n
    if isnumeric(list{i})
      list{i} = sprintf(lParams.doubleFormat, list{i});
    elseif isinteger(list{i})
      list{i} = sprintf(lParams.integerFormat, list{i});
    endif
  endfor

  S = strcat(strjoin(list(1:n), delimiter), lPrefix);
endfunction
