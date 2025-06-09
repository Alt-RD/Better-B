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
%  Copyright (c) 2022 Montpellier-University
%  Copyright (c) 2022-2025 AltRD: Emmanuel Ruffio
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
function v = HT_Material_GetLambda(mat, dir)
  assert(HT_CheckType(mat, "material"), 'Invalid material object');
  assert(isfield(mat, 'lambda'));
  if nargin < 2, dir = []; endif

  if isempty(dir)
    v = NA(3, numel(mat));

    for i=1:numel(mat)
      if isscalar(mat(i).lambda)
        v(:,i) = repmat(mat(i).lambda, 3, 1);
      else
        v(:,i) = mat(i).lambda;
      endif
    endfor
  else
    assert((dir >= 1) && (dir <= 3), 'Invalid direction');
##    assert(dir <= numel(mat.lambda), 'Invalid material conductivity');
##    v = mat.lambda(dir);

    v = NA(numel(mat), 1);
    for i=1:numel(mat)
      if isscalar(mat(i).lambda)
        v(i) = mat(i).lambda;
      else
        v(i) = mat(i).lambda(dir);
      endif
    endfor
  endif
endfunction
