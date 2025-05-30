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
% Voir Radiation View factor (p.28)
% Return view factor F12
% f1: Positions du rectangle 1
%   > f1 = xm xp ym yp
% f2: Positions du rectangle 2
%   > f2 = xm xp ym yp
% xm and xp are the distance to the common axis.
% ym and yp are the position along the common axis.
function F = HT_ViewFactor_Rect2RectAngle(f1, f2)
  ind1 = [repmat([1 2], 8, 1)(:)  repmat([3 4 3 4], 4, 1)(:)];
  ind2 = [repmat([1;1;2;2], 4, 1) repmat([3;4], 8, 1)];
  K = [1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1];
  F = zeros(size(f1,1), 1);

  s = warning("query", "Octave:divide-by-zero");
  warning("off", "Octave:divide-by-zero");

  for i=1:16
    x1 = f1(:, ind1(i, 1));
    y1 = f1(:, ind1(i, 2));
    x2 = f2(:, ind2(i, 1));
    y2 = f2(:, ind2(i, 2));

    Csqr = max(x1.^2+x2.^2, 0); % Avoid small negative value
    C = sqrt(Csqr);
    Dy = (y1-y2);
    DySqr = Dy.^2;

    A = atan(Dy ./ C);  %    A(D > 1000) = pi/2 - atan(1./D);
    A(C == 0) = pi/2;

    L = (Csqr-DySqr).*log(Csqr+DySqr);
    L( (Csqr+DySqr)==0.0 ) = 0.0;

    B = Dy .* C .* A - ...
         0.25*L;

    F += K(i) * B;
  endfor

  A = (f1(:,2)-f1(:,1)) .* (f1(:,4)-f1(:,3));
  F ./= 2*pi*A;

  warning(s);
endfunction
