%  This file is part of project HiveModel.
%  This work was supported by the Better-B project, which has received funding
%  from the European Union, the Swiss State Secretariat for Education, Research
%  and Innovation (SERI) and UK Research and Innovation (UKRI) under the UK
%  government's Horizon Europe funding guarantee (grant number 10068544). Views
%  and opinions expressed are however those of the author(s) only and do not
%  necessarily reflect those of the European Union, European Research Executive
%  Agency (REA), SERI or UKRI. Neither the European Union nor the granting
%  authorities can be held responsible for them.
%
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2023-2025: CoActions-AltRD-Emmanuel Ruffio
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

lUserPlotOptions.currentFigureIndex = 1;


% Plot the sun position with time if requested
if isfield(lUserPlotOptions, 'plotSunPosition') && lUserPlotOptions.plotSunPosition
  lFigureIndex = lUserPlotOptions.figureIndex(lUserPlotOptions.currentFigureIndex);
  lUserPlotOptions.currentFigureIndex++;
  figure(lFigureIndex);
  clear lFigureIndex;

  clf;
  plot(timeVector, 180*(pi/2-zenithAngle)/pi, 'displayname', 'Elevation Angle'); hold on;
  plot(timeVector, 180*azimuthAngle/pi, 'displayname', 'Azimuth Angle');
  xlabel('Time', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  ylabel('Angle', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  set(gca, 'fontsize', lPlotOptions.axisFontSize, 'fontweight', lPlotOptions.axisFontWeight);
  title('Sun position vs time', 'fontsize', lPlotOptions.titleFontSize, 'fontweight', lPlotOptions.titleFontWeight);
  hl = legend();

##  plot3(lSunDirectionMatrix(1,:), lSunDirectionMatrix(2,:), lSunDirectionMatrix(3,:));
  grid on;
endif

% Plot Cos of angle between sides and sun
if isfield(lUserPlotOptions, 'plotSunSideDirection') && lUserPlotOptions.plotSunSideDirection
  lFigureIndex = lUserPlotOptions.figureIndex(lUserPlotOptions.currentFigureIndex);
  lUserPlotOptions.currentFigureIndex++;
  figure(lFigureIndex);
  clear lFigureIndex;

  clf;
  subplot(1,2,1);
  h = plot(timeVector/3600, cosBeta); hold on;
  lDisplayNames = {'0deg (north)', '45deg', '90deg (west)', '135deg', '180deg (south)', '225deg', '270deg (east)', '315deg', 'Top (roof)'};
  lDisplayStyle = {'-', '-', '-', '-', ':', ':', ':', ':', '-.'};
  arrayfun(@(i) set(h(i), 'displayname', lDisplayNames{i}, 'linestyle', lDisplayStyle{i}, 'linewidth', 2), 1:columns(cosBeta));
  xlabel('Time [hour]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  ylabel('Angle', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  set(gca, 'fontsize', lPlotOptions.axisFontSize, 'fontweight', lPlotOptions.axisFontWeight);
  title('Cosine of angle between sides and sun', 'fontsize', lPlotOptions.titleFontSize, 'fontweight', lPlotOptions.titleFontWeight);
  hl = legend();
  grid on;

  subplot(1,2,2);
  h = plot(timeVector/3600, lSunIrradianceMatrix); hold on;
  lDisplayNames = {'0deg (north)', '45deg', '90deg (west)', '135deg', '180deg (south)', '225deg', '270deg (east)', '315deg', 'Top (roof)'};
  lDisplayStyle = {'-', '-', '-', '-', ':', ':', ':', ':', '-.'};
  arrayfun(@(i) set(h(i), 'displayname', lDisplayNames{i}, 'linestyle', lDisplayStyle{i}, 'linewidth', 2), 1:columns(cosBeta));
  xlabel('Time [hour]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  ylabel('Irradiance [W/m²]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
  set(gca, 'fontsize', lPlotOptions.axisFontSize, 'fontweight', lPlotOptions.axisFontWeight);
  title('Irradiance on each sides', 'fontsize', lPlotOptions.titleFontSize, 'fontweight', lPlotOptions.titleFontWeight);
  hl = legend();
  set(hl, 'location', 'northwest');
  grid on;

endif


