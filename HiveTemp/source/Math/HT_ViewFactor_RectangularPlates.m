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
% Return the view factor of 2 rectangular plates of size W1, W2
% separated by distance H (cf p.23, Radiation View Factors)
% which are perfectly aligned (face to face, not coplanar)
function F = HT_ViewFactor_RectangularPlates(W1, W2, H)
  x = W1./H;
  y = W2./H;
  x1Sqr = 1+x.^2;
  y1Sqr = 1+y.^2;
  x1 = sqrt(x1Sqr);
  y1 = sqrt(y1Sqr);

  F = log(x1Sqr*y1Sqr ./ (x1Sqr+y1Sqr-1)) + ...
      2*x.*(y1.*atan(x./y1) - atan(x)) + ...
      2*y.*(x1.*atan(y./x1) - atan(y));
  F ./= (pi*x.*y);
endfunction
