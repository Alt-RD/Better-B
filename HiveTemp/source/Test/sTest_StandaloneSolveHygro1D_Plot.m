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
%  Copyright (c) 2023-2025 AltRD-Emmanuel Ruffio
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

assert(exist('WMatNL', 'var') == 1, 'Script sTest_StandaloneSolveHygro1D must be executed first');

figure(1);
clf;
tPlot = (0.0:0.2:1)*lParamsNLCase.t(end);
x = lInfosNL.nodePosition * lParamsNLCase.length * 100;
xcenter = mean(x);
lSecondToDay = @(s) s/3600/24;

hList = [];

for k=1:numel(tPlot)
  ind = find(tPlot(k) >= lParamsNLCase.t, 1, 'last');
  y = WMatNL(:,ind);
  h = plot(x, y,  ...
    'linewidth', 2, ...
    'displayname', sprintf('t=%.1f', lSecondToDay(lParamsNLCase.t(ind)))); hold on;

  text(xcenter, interp1(x,y,xcenter)+0.0002, sprintf('t=%.0f days', lSecondToDay(lParamsNLCase.t(ind))), ...
    'fontsize', 14, 'fontweight', 'bold');

  hList = [hList h];
endfor

for k=1:numel(tPlot)
  ind = find(tPlot(k) >= lParamsNLCase.t, 1, 'last');
  y = WMatL(:,ind);
  plot(x, y,  ...
    'linewidth', 2, ...
    'linestyle', '--', ...
    'color', get(hList(k), 'color'), ...
    'displayname', sprintf('t=%.1f', lSecondToDay(lParamsNLCase.t(ind)))); hold on;

  hList = [hList h];
endfor

grid on;
xlabel('Position [cm]', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Water content profiles [-]', 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'fontsize', 16, 'fontweight', 'bold');

tmp = [lParamsNLCase.winit lWeqVector];
ylim([min(tmp) max(tmp)]);

figure(2);
clf;
plot(lParamsNLCase.t/3600, WMatNL(1,:), 'linewidth', 2.0); hold on;
plot(lParamsNLCase.t/3600, WMatNL(end,:), '+'); hold on;
plot(lParamsNLCase.t/3600, WMatNL(1+(lParamsNLCase.n-1)/2,:), '-k', 'displayname', 'center', 'linewidth', 2.0); hold on;
plot(lParamsNLCase.t/3600, WMatL(1+(lParamsNLCase.n-1)/2,:), '--k', 'displayname', 'center', 'linewidth', 2.0); hold on;
grid on;
xlabel('Time [hour]', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Water content [-]', 'fontsize', 16, 'fontweight', 'bold');
title('Water content at the center', 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'fontsize', 16, 'fontweight', 'bold');



