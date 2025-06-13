%  This file is part of project HiveTemp.
%
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

assert(isvarname('T') && isvarname('U') && isvarname('Tnodes'), 'Script sHiveModel_Datant.m must be started first');

% ======================= Display temperature evolution ========================
if lPlotOptions.onlyOneFigure, clf; subplot(2,3,[2 3]); else figure(2); clf; endif
% Extract the temperature evolution of node called "nbody"
tplot = 0:600:t(end);
lPlotNodesInfos = struct( 'nodeName', { 'nAirExt', ...
                                        'nHBBot', ...
                                        lMod_HiveBodyVolumes(1).nodes{:}, ...
                                        lMod_HiveBodyVolumes(2).nodes{:}
                                        }, ...
                          'displayName', {'Air ext.', ...
                                          'Hive body bottom', ...
                                          lMod_HiveBodyVolumes(1).nodes{:}, ...
                                          lMod_HiveBodyVolumes(2).nodes{:}}, ...
                          'lineStyle', {'-', ...
                                        '-', ...
                                        repmat({'--'}, 1, numel(lMod_HiveBodyVolumes(1).nodes)){:}, ...
                                        repmat({'--'}, 1, numel(lMod_HiveBodyVolumes(2).nodes)){:}}, ...
                          'lineColor', {'black', ...
                                        'red', ...
                                        repmat({'blue'}, 1, numel(lMod_HiveBodyVolumes(2).nodes)){:}, ...
                                        repmat({'blue'}, 1, numel(lMod_HiveBodyVolumes(2).nodes)){:}}, ...
                          'lineWidth', 4);

lPlotVolumeInfos = struct('volumes', { lMod_HiveBodyVolumes }, ...
                     'displayName', 'Hive body air', ...
                    'lineStyle', '-', ...
                    'lineColor', {'blue'}, ...
                    'lineWidth', 4);

% Build the list of faces whose average temperature must be plotted
lPlotFacesList = { HT_Face_Find(lMod_SideWall1_facesCell{1}, 'orientation', lHiveParams.globalAxis(:,X), 'numberCheck', 1, 'returnType', true), ...
                   HT_Face_Find(lMod_SideWall1_facesCell{1}, 'orientation', -lHiveParams.globalAxis(:,X), 'numberCheck', 1, 'returnType', true) };
lPlotFacesInfos = struct( 'faces', lPlotFacesList, ...
                          'displayName', {'SW1_XM_internal', 'SW1_XM_external'}, ...
                          'lineStyle', '-', ...
                          'lineColor', {'blue', 'green'}, ...
                          'lineWidth', 4);

lNodeTemp = HT_Result_GetNodeTemperature(...
                    {lPlotNodesInfos.nodeName},   ... List of nodes
                    T,                      ... Temperature vector
                    Tnodes,                 ... Model modes list
                    'time', tplot,          ... Display time vector
                    'timevector', t);       ... Time Vector

lFaceTemp = HT_Result_GetFaceTemperature(...
                    {lPlotFacesInfos.faces},   ... List of nodes
                    T,                      ... Temperature vector
                    Tnodes,                 ... Model modes list
                    'time', tplot,          ... Display time vector
                    'timevector', t, ...    ... Time Vector
                    'operation', 'array');

lVolumeTemp = HT_Result_GetVolumeTemperature(...
                    lPlotVolumeInfos.volumes, ...
                    T,                      ... Temperature vector
                    Tnodes,                 ... Model modes list
                    'time', tplot,          ... Display time vector
                    'timevector', t, ...    ... Time Vector
                    'operation', 'average');

##Test = HT_Result_GetNodeTemperature(...
##                    {lMod_HiveBodyVolumes(2).nodes},   ... List of nodes
##                    T,                      ... Temperature vector
##                    Tnodes,                 ... Model modes list
##                    'time', tplot,          ... Display time vector
##                    'timevector', t);       ... Time Vector

for i=1:numel(lPlotNodesInfos)
  plot(tplot/3600,  lNodeTemp(i,:),               ... Temperature data
                    lPlotNodesInfos(i).lineStyle, ... Line style
                    'linewidth', lPlotNodesInfos(i).lineWidth, ... Line width
                    'color', lPlotNodesInfos(i).lineColor,     ... Line color
                    'displayname', lPlotNodesInfos(i).displayName); hold on;
endfor

for i=1:numel(lPlotFacesInfos)
  plot(tplot/3600,  lFaceTemp(i,:),               ... Temperature data
                    lPlotFacesInfos(i).lineStyle, ... Line style
                    'linewidth', lPlotFacesInfos(i).lineWidth, ... Line width
                    'color', lPlotFacesInfos(i).lineColor,     ... Line color
                    'displayname', lPlotFacesInfos(i).displayName); hold on;
endfor

for i=1:numel(lPlotVolumeInfos)
  plot(tplot/3600,  lVolumeTemp(i,:),               ... Temperature data
                    lPlotVolumeInfos(i).lineStyle, ... Line style
                    'linewidth', lPlotVolumeInfos(i).lineWidth, ... Line width
                    'color', lPlotVolumeInfos(i).lineColor,     ... Line color
                    'displayname', lPlotVolumeInfos(i).displayName); hold on;
endfor

set(gca, 'fontsize', lPlotOptions.axisFontSize, 'fontweight', lPlotOptions.axisFontWeight);
ylabel('Temperature [degC]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
xlabel('Time [hour]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
title('Temperature of hive body and of int./ext. face of each side', ...
    'fontsize', lPlotOptions.titleFontSize, 'fontweight', lPlotOptions.titleFontWeight);
hl = legend();
set(hl, 'fontsize', 10);
grid on;

##% ================ Display temperature profils inside Xp wall ==================
##if lPlotOptions.onlyOneFigure, subplot(2,3,4); else figure(3); clf; endif
##% Extract temperatures for each nodes of XP wall. The list of nodes was saved
##% in <lWallNodes{FaceXP}>. The time is interpolated and is regularly spaced.
####tprofil = logspace(0, log10(t(end)), 10); %0:600:t(end); % [s]
##tprofil = 0:7200:t(end); % [s]
##nNodes = numel(lWallNodes{FaceXP});
##Twallprofil = HT_Result_GetNodeTemperature(lWallNodes{FaceXP}, T, Tnodes, ...
##                                          'time', tprofil, ... Time of different profils
##                                          'timevector', t);
##for k=1:size(Twallprofil,2)
##  h = plot(Twallprofil(:,k), 'linewidth', 2.0); hold on;
##  disp(get(h, 'color'));
##  text(nNodes/2, Twallprofil(idivide(nNodes,2),k), ...
##        sprintf('t=%.1fh', tprofil(k)/3600),...
##        'fontsize', 14, ...
##        'fontweight', 'bold');
##endfor
##
##set(gca, 'fontsize', 14);
##xlabel('Position [node]');
##ylabel('Temperature [degC]');
##title('Temperature profils inside hive wall XP (north)');
##grid on;
##
##% ================ Display temperature profils inside Xm wall ==================
##if lOnlyOneFigure, subplot(2,3,5); else figure(4); clf; endif
##% Extract temperatures for each nodes of XP wall. The list of nodes was saved
##% in <lWallNodes{FaceXP}>. The time is interpolated and is regularly spaced.
####tprofil = logspace(0, log10(t(end)), 10); %0:600:t(end); % [s]
##tprofil = 0:7200:t(end); % [s]
##nNodes = numel(lWallNodes{FaceXM});
##Twallprofil = HT_Result_GetNodeTemperature(lWallNodes{FaceXM}, T, Tnodes, ...
##                                          'time', tprofil, ... Time of different profils
##                                          'timevector', t);
##for k=1:size(Twallprofil,2)
##  plot(Twallprofil(:,k), 'linewidth', 2.0); hold on;
##  text(nNodes/2, Twallprofil(idivide(nNodes,2),k), ...
##        sprintf('t=%.1fh', tprofil(k)/3600),...
##        'fontsize', 14, ...
##        'fontweight', 'bold');
##endfor
##
##set(gca, 'fontsize', 14);
##xlabel('Position [node]');
##ylabel('Temperature [degC]');
##title('Temperature profils inside hive wall XM (south)');
##grid on;

% =================== Display solar irradiance for each face ===================
if lPlotOptions.onlyOneFigure, subplot(2,3,6); else figure(5); clf; endif
hold on;
lFaceName = { 'XM', 'XP', 'YM', 'YP', 'ZM', 'ZP' };
for i=1:6
  h = plot(t/3600,  lHiveBodySWExpositionMat(:,i),              ...
                    'linewidth', 3.0,                 ...
                    'color', lPlotColors.faces(i,:),  ...
                    'displayname', sprintf('Sidewalls %s', lFaceName{i}));
endfor
for i=1:6
  h = plot(t/3600,  lRoofExpositionMat(:,i),              ...
                    'linewidth', 3.0,                 ...
                    'color', lPlotColors.faces(i,:),  ...
                    'linestyle', ':', ...
                    'displayname', sprintf('Roof %s', lFaceName{i}));
endfor
##plot(t/3600, lSolarIrradiance, '--k', 'linewidth', 3.0, 'displayname', 'Solar irradiance');
set(gca, 'fontsize', lPlotOptions.axisFontSize, 'fontweight', lPlotOptions.axisFontWeight);
ylabel('Solar irradiance [W/m²]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
xlabel('Time [hour]', 'fontsize', lPlotOptions.labelFontSize, 'fontweight', lPlotOptions.labelFontWeight);
legend();
title('Direct solar irradiance incoming on each hive faces (wood part)', ...
    'fontsize', lPlotOptions.titleFontSize, 'fontweight', lPlotOptions.titleFontWeight);
grid on;

