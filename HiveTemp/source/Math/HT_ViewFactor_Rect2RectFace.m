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
%
% Voir Radiation View factor (p.23)
% Between non coplanar rectangle separated by a distance <z>, but not necessarily
% aligned (contrary to HT_ViewFactor_RectangularPlates)
% Return view factor F12
% f1: Positions du rectangle 1
%   > f1 = xm xp ym yp
% f2: Positions du rectangle 2
%   > f2 = xm xp ym yp
% Distance verticale entre les faces
function F = HT_ViewFactor_Rect2RectFace(f1, f2, z)
  ind1 = [repmat([1 2], 8, 1)(:)  repmat([3 4 3 4], 4, 1)(:)];
  ind2 = [repmat([1;1;2;2], 4, 1) repmat([3;4], 8, 1)];
  K = [1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1];
  F = zeros(size(f1,1), 1);

  % Every faces (row of f1 and f2) are computed simultaneously
  for i=1:16
    x1 = f1(:, ind1(i, 1));
    y1 = f1(:, ind1(i, 2));
    x2 = f2(:, ind2(i, 1));
    y2 = f2(:, ind2(i, 2));

    u = x1-x2;
    uSqr = u.^2;
    v = y1-y2;
    vSqr = v.^2;
    p = sqrt(uSqr + z^2);
    q = sqrt(vSqr + z^2);

    F += K(i) * (v.*p.*atan(v./p) + ...
                 u.*q.*atan(u./q) - ...
                 z^2/2*log(uSqr+vSqr+z^2));
  endfor

  A = (f1(:, 2)-f1(:, 1)).*(f1(:, 4)-f1(:, 3)); % Surface area of every faces
  F ./= 2*pi*A;
endfunction
