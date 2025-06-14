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
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2025 AltRD-Emmanuel Ruffio
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
if lPlotOptions.onlyOneFigure, clf; else figure(2); clf; endif

for i=1:nParam
  if lPlotOptions.onlyOneFigure, clf; subplot(1,nParam,i); else figure(1+i); clf; endif

  disp(sprintf('Plotting sensitivity of parameter <%s>', XInfos(i).name));

  % Extract the temperature evolution of node called "nbody"
  tplot = 0:600:t(end);
  lPlotNodesInfos = struct( 'nodeName', {'nAirExt', 'nHBBot',}, ...
                            'displayName', {'Air ext.', 'Hive body bottom'}, ...
                            'lineStyle', '-', ...
                            'lineColor', {'black', 'red'}, ...
                            'lineWidth', 4);

  % Build the list of faces whose average temperature must be plotted
  lPlotFacesList = { HT_Face_Find(XInfos(i).faces.sideWall1{1}, 'orientation', lHiveParams.globalAxis(:,X), 'numberCheck', 1, 'returnType', true), ...
                     HT_Face_Find(XInfos(i).faces.sideWall1{1}, 'orientation', -lHiveParams.globalAxis(:,X), 'numberCheck', 1, 'returnType', true) };
  lPlotFacesInfos = struct( 'faces', lPlotFacesList, ...
                            'displayName', {'SW1_XM_internal', 'SW1_XM_external'}, ...
                            'lineStyle', '-', ...
                            'lineColor', {'blue', 'green'}, ...
                            'lineWidth', 4);

  lPlotVolumeInfos = struct('volumes', { XInfos(i).volumes.hiveBody }, ...
                       'displayName', 'Hive body air', ...
                      'lineStyle', '-', ...
                      'lineColor', {'blue'}, ...
                      'lineWidth', 4);

  lNodeTemp = HT_Result_GetNodeTemperature(...
                      {lPlotNodesInfos.nodeName},   ... List of nodes
                      XInfos(i).X,                  ... Sensitivity matrix
                      XInfos(i).nodes,              ... Model modes list
                      'time', tplot,                ... Display time vector
                      'timevector', t);             ... Time Vector

  lFaceTemp = HT_Result_GetFaceTemperature(...
                      {lPlotFacesInfos.faces},      ... List of nodes
                      XInfos(i).X,                  ... Temperature vector
                      XInfos(i).nodes,              ... Model modes list
                      'time', tplot,                ... Display time vector
                      'timevector', t,              ... Time Vector
                      'operation', 'array');

  lVolumeTemp = HT_Result_GetVolumeTemperature(...
                      lPlotVolumeInfos.volumes, ...
                      XInfos(i).X,                    ... Temperature vector
                      XInfos(i).nodes,                ... Model modes list
                      'time', tplot,                  ... Display time vector
                      'timevector', t,                ... Time Vector
                      'operation', 'average');

  for k=1:numel(lPlotNodesInfos)
    plot(tplot/3600,  lNodeTemp(k,:),               ... Temperature data
                      lPlotNodesInfos(k).lineStyle, ... Line style
                      'linewidth', lPlotNodesInfos(k).lineWidth, ... Line width
                      'color', lPlotNodesInfos(k).lineColor,     ... Line color
                      'displayname', lPlotNodesInfos(k).displayName); hold on;
  endfor

  for k=1:numel(lPlotFacesInfos)
    plot(tplot/3600,  lFaceTemp(k,:),               ... Temperature data
                      lPlotFacesInfos(k).lineStyle, ... Line style
                      'linewidth', lPlotFacesInfos(k).lineWidth, ... Line width
                      'color', lPlotFacesInfos(k).lineColor,     ... Line color
                      'displayname', lPlotFacesInfos(k).displayName); hold on;
  endfor


  set(gca, 'fontsize', 14);
  xlabel('Temps [hour]');
  ylabel('Temperature [degC]');
  title('Temperature of hive body and of int./ext. face of each side');
  hl = legend();
  set(hl, 'fontsize', 10);
  grid on;

  title(sprintf('Sensitivity <%s>', XInfos(i).name));


endfor

##lPlotInfos = {'nHBBot'      , 'Hive body (bot)' , '-' , 'red'                   , 2; ... List of nodes whose temperature
##              'nHBTop'      , 'Hive body (top)' , '--', 'blue'                  , 2};
##
##if lPlotConfig.onlyOneFigure
##  figure(1);
##  clf;
##endif
##
##lTimeScale = [1 1 60 60 3600];
##lTimeShortName = {'sec', 'sec', 'min', 'min' 'hour'};
##lTimeUnitIndex = find(strcmpi(lPlotConfig.timeUnit, {'sec', 'second', 'min', 'minute', 'hour'}));
##tplot = (t(1):lPlotConfig.fixedTimeStep:t(end));
##
##nInfos = rows(lPlotInfos);
##
##for i=1:nPlotParam
##  if ~lPlotConfig.onlyOneFigure
##    figure(i);
##    clf;
##  endif
##
##  ind = lPlotParamIndex(i);
##  lPlotConf = lPlotParamList(ind,:);
##  lPlotName = lPlotConf{2};
##  lPlotSettings = lPlotConf{3};
##
##  lPlotTemp = HT_Result_GetNodeTemperature(...
##                      lPlotInfos(:,1),          ... List of nodes
##                      lSensitivityMatrix{ind},    ... Temperature vector
##                      Tnodes,                   ... Model modes list
##                      'time', tplot,            ... Display time vector
##                      'timevector', t);         ... Time Vector
##
##  h = plot(tplot/ lTimeScale(lTimeUnitIndex), lPlotTemp);
##
##  for k=1:nInfos
##    set(h(k),    'linestyle', lPlotInfos{k,3}, ...
##                 'color', lPlotInfos{k,4}, ...
##                 'linewidth', lPlotInfos{k,5}, ...
##                 'displayname', sprintf('%s.%s', lPlotName, lPlotInfos{k,2}));
##  endfor
##
##  if ~isempty(lPlotSettings), set(h, lPlotSettings{:}); endif
##
####  for i=1:idivide(numel(lPlotSettings), int32(2))
####    set(h, lPlotSettings{2*i-1}, lPlotSettings{2*i});
####  endfor
##  if ~lPlotConfig.onlyOneFigure
##    xlabel(sprintf('Time [%s]', lTimeShortName{lTimeUnitIndex}), 'fontsize', 20, 'fontweight', 'bold');
##    ylabel('Sensitivity [K]', 'fontsize', 20, 'fontweight', 'bold');
##    hl = legend();
##    set(hl, 'fontsize', 14, 'fontweight', 'bold');
##    grid on;
##  endif
##
##endfor
##
##if lPlotConfig.onlyOneFigure
##  xlabel('Time [min]', 'fontsize', 20, 'fontweight', 'bold');
##  ylabel('Sensitivity [K]', 'fontsize', 20, 'fontweight', 'bold');
##  hl = legend();
##  set(hl, 'fontsize', 14, 'fontweight', 'bold');
##  grid on;
##endif
##

