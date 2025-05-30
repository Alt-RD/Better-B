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
% Return the view factor F12 of 2 rectangular plates adjacent
% of size HxL and WxL (cf p.27, Radiation View Factors)
% Warning the face 1 (WxL) is the horizontal rectangle and face 2
% (HxL) the vertical one
function F = HT_ViewFactor_RectangularPlatesAngle(H, W, L)
  if numel(H) == 1 && numel(W) == 1
    h = H/L;
    w = W/L;
    hSqr = h.^2;
    wSqr = w.^2;

    t = (1 + hSqr + wSqr) ./ (hSqr + wSqr);
    b = t .* wSqr ./ (1 + wSqr);
    c = t .* hSqr ./ (1 + hSqr);
    a = (1 + hSqr) .* (1 + wSqr) ./ (1+hSqr+wSqr);

    F = h .* atan(1./h) + w .* atan(1./w) - ...
        sqrt(hSqr + wSqr) .* atan(1 ./ sqrt(hSqr + wSqr)) + ...
        0.25 * log(a.*b.^(wSqr).*c.^(hSqr));
    F ./= (pi*w);
  elseif numel(H) == 2 && numel(W) == 1
    % The second face is not connected to face 1
    assert(H(2) > H(1));
    F = HT_ViewFactor_RectangularPlatesAngle(H(2), W, L) - ...
        HT_ViewFactor_RectangularPlatesAngle(H(1), W, L);
  elseif numel(H) == 1 && numel(W) == 2
    assert(W(2) > W(1));
    F = (H/(W(2)-W(1)))*( HT_ViewFactor_RectangularPlatesAngle(W(2), H, L) - ...
                          HT_ViewFactor_RectangularPlatesAngle(W(1), H, L));
  elseif numel(H) == 2 && numel(W) == 2
    A22p = H(2);

    F = H(2)/(W(2)-W(1)) * HT_ViewFactor_RectangularPlatesAngle(W, H(2), L) - ...
        H(1)/(W(2)-W(1)) * HT_ViewFactor_RectangularPlatesAngle(W, H(1), L);
  else
    error('invalid parameter');
  endif
endfunction
