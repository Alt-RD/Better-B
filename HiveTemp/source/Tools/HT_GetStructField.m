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
% Returns the field <fieldName> of every object x
function V = HT_GetStructField(x, fieldName, reduce)
  if nargin < 3, reduce = true; endif;
  assert(islogical(reduce));

  if iscell(x)
    assert(all(cellfun(@(v) isstruct(v), x)), 'x must be a structure or a cell array of structure');
    V = cellfun(@(v) HT_GetStructField(v, fieldName), x, 'uniformoutput', false);

    if reduce, V = cell2mat(V); endif;
  elseif isstruct(x)
    V = arrayfun(@(v) getfield(v, fieldName), x, 'Uniformoutput', false);
  else
    error('HT_GetStructField::Invalid object');
  endif
endfunction
