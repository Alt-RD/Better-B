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
function V = HT_AddToTimeVec(V, N, unit)
  assert(any(strcmpi(unit, {'year', 'month', 'day', 'hour', 'minute', 'second'})));
  assert((numel(V) == 6) || (numel(V) == 3));
  assert(numel(N) == 1);

  dayPerMonthf = @(year) [31 28+is_leap_year(year) 31 30 31 30 31 31 30 31 30 31];

  if strcmpi(unit, 'year')
    V(1) += N;
  elseif strcmpi(unit, 'month')
    V(2) += N;
    lAddYear = floor(V(2)-1)/12;
    V(2) -= lAddYear * 12;
    V(1) += lAddYear;
  elseif strcmpi(unit, 'day')
    V(3) += N;

    while (V(3) > dayPerMonthf(V(1))(V(2)))
      V(3) -= dayPerMonthf(V(1))(V(2));
      V(2) += 1;
      if V(2) > 12
        V(2) -= 12;
        V(1)++;
      endif
    endwhile
  elseif strcmpi(unit, 'hour')
    V(4) += N;
    lAddDay = floor(V(4)/24);
    V = HT_AddToTimeVec(V, lAddDay, 'day');
    V(4) -= lAddDay * 24;
  else
    error('Not implemented');
  endif
endfunction
%!test
%! V = [2023 02 16 16 17 00];
%! assert(HT_AddToTimeVec(V, 1, 'day'), [2023 02 17 16 17 00]);
%!
%!test
%! V = [2023 02 16 16 17 00];
%! assert(HT_AddToTimeVec(V, 13, 'day'), [2023 03 1 16 17 00]);
%!
%!test
%! V = [2023 02 16 16 17 00];
%! assert(HT_AddToTimeVec(V, 44, 'day'), [2023 04 1 16 17 00]);
%!
%!test
%! V = [2023 02 16 16 17 00];
%! assert(HT_AddToTimeVec(V, 11, 'month'), [2024 01 16 16 17 00]);

