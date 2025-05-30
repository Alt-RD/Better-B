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
% Returns an array by merging A and B
% Y(i) = A(i) if not NAN, else B(i)
% Error if A(i) and B(i) not NaN and different
function R = HT_MergeArray(A, B, priority)
  if nargin < 3, priority = 0; endif

  R = NA(size(A));

  assert(numel(A) == numel(B));

  lValidA = ~(isnan(A) | isna(A));
  lValidB = ~(isnan(B) | isna(B));

  if priority == 0
    assert(~any(lValidA & lValidB), 'Some elements of the array can not be merged.');
  elseif priority == 1
    lValidB(lValidA) = false;
  elseif priority == 2
    lValidA(lValidB) = false;
  else
    error('Invalid value of priority');
  endif

  R(lValidA) = A(lValidA);
  R(lValidB) = B(lValidB);
endfunction
