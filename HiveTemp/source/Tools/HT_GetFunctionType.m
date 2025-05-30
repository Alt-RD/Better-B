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

% Check the type of a function handle
% time varying or non linear time varying
% Returns the number of parameter for each function
function T = HT_GetFunctionType(funcList)
  if ~iscell(funcList), funcList = { funcList }; endif;

  nf = numel(funcList);
  T = zeros(nf, 1);

  for i=1:nf
    assert(is_function_handle(funcList{i}), 'Handle is not a function');

    funStr = functions(funcList{i}).function;
    ind1 = strfind(funStr, '(')(1);
    ind2 = strfind(funStr, ')')(1);
    paramCount = numel(strsplit(funStr((ind1+1):(ind2-1)), ','));

    assert(any(paramCount == [2 3]), 'Functions must have 2 (t,i) or 3 (T,t,i) parameters');
    T(i) = paramCount;
  endfor
endfunction
