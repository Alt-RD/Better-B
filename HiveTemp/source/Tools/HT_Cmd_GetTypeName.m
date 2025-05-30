%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022: Montpellier University / CoActions-AltRD-Emmanuel Ruffio
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
function S = HT_Cmd_GetTypeName(type)
  HT_ImportConstants();

  switch (type)
    case HT_BTYPE_T_VARIABLE_UNIFORM
      S = '(T) Variable and uniform';
    case HT_BTYPE_T_CONSTANT_UNIFORM
      S = '(T) Constant and uniform';
    case HT_BTYPE_T_INTERPOLATION_ARRAY
      S = '(T) Interpolation array';

    case HT_BTYPE_FLUX_VARIABLE_UNIFORM
      S = '(Flux) Variable and uniform';
    case HT_BTYPE_FLUX_VARIABLE_NON_UNIFORM
      S = '(Flux) Variable and non uniform';
    case HT_BTYPE_FLUX_CONSTANT_UNIFORM
      S = '(Flux) Constant and uniform';
    case HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM
      S = '(Flux) Constant and non uniform';
    case HT_BTYPE_FLUX_INTERPOLATION_ARRAY
      S = '(Flux) Interpolation array';

    otherwise
      S = 'unknown';
  endswitch

endfunction
