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
% This program computes the water vapour transfer through the hive walls.
% To get more consistent values, the hygro-thermal climate inside the hive is
% given by hive measurements. The same computations could be done by building a
% model of the hive but it would involves parameters that are difficult to estimate.
%
% THe water transfer through the hive walls is dependant on the temperature of
% both faces of the hive walls. As a result, a 1D thermal model is defined.
clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../../HiveTemp/source/';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));
% Init the library
HT_Init();

% Source files
lFilePathEnv = make_absolute_filename('./Input/2024_12_23/2025_01_02_Env.csv');
lFilePath = make_absolute_filename('./Input/2024_12_23/2025_01_02_TestHive.csv');
lLine = '';

[lPath lName] = fileparts(lFilePath);
lCacheFile = strcat('./Cache/',lName, '.mat');

[lPathEnv lNameEnv] = fileparts(lFilePathEnv);
lCacheFileEnv = strcat('./Cache/',lNameEnv, '.mat');

% Compute the average water content in the hive body based on humidity and
% temperature sensors
% Define constant for convenience
EPOCH = 1;
TOP_RH = 2;
TOP_T = 3;
BOT_RH = 4;
BOT_T = 5;
SOUTH_EXT_T = 6;
SOUTH_EXT_PTAT = 7;
EXT_T = 8;
EXT_RH = 9;
EXT_T2 = 10;

lSensors = struct('name', { 'Rtc.epoch [day]', ...
                            'Test.body.top.HR [%]', ...
                            'Test.body.top.T [deg]', ...
                            'Test.body.bottom.HR [%]', ...
                            'Test.body.bottom.T [deg]', ...
                            'Test.camera.ext.south.T [deg]', ...
                            'Test.camera.ext.south.PTAT.T [deg]' }, ...
                  'data', [], ...
                  'scale', { 1, 0.01, 1, 0.01, 1, 1, 1 },
                  'offset', { 0, 0, 273.15, 0, 273.15, 273.15, 273.15}, ... Convert deg C to Kelvin C
                  'color', {[], ...
                            [224 137 139]/255,...
                            [130 151 238]/255, ...
                            [170 45 48]/255, ...
                            [13 35 126]/255, ...
                            [255 127 39]/255, ...
                            [0 255 255]/255 });

lSensorsEnv = struct('name', {  'Air1.T [deg]', ...
                                'Air1.HR [%]', ...
                                'Air2.T [deg]'}, ...
                  'data', [], ...
                  'scale', { 1, 0.01, 1 },
                  'offset', { 273.15, 0, 273.15}, ... Convert deg C to Kelvin C
                  'color', {[0 0 0]/255,...
                            [ 0 0 0]/255, ...
                            [0 0 0]/255});

[D] = HT_ReadCsvFile(lFilePath,         'cacheFile', lCacheFile, ...
                                        'columns', {lSensors.name}, ...
                                        'forceReload', false, ...
                                        'subSampling', 3);

[DEnv] = HT_ReadCsvFile(lFilePathEnv,   'cacheFile', lCacheFileEnv, ...
                                        'columns', {lSensorsEnv.name}, ...
                                        'forceReload', false, ...
                                        'subSampling', 3);

lSubSampling = 3;
% Load the data from matrix <D> and scale the data with appropriate coefficient
% if necessary (ex: HR % -> HR)
for i=1:numel(lSensors)
  lSensors(i).data = lSensors(i).scale * D{1,i};
  lSensors(i).data += lSensors(i).offset;
endfor

for i=1:numel(lSensorsEnv)
  lSensorsEnv(i).data = lSensorsEnv(i).scale * DEnv{1,i};
  lSensorsEnv(i).data += lSensorsEnv(i).offset;
endfor

lSensors = [lSensors, lSensorsEnv];
clear lSensorEnv;

% Convert time unit from day to hour and set the initial time
lTimeVec = lSensors(EPOCH).data - floor(lSensors(EPOCH).data(1));
lTimeVec *= 24;

% Compute the vapour partial pressure at the top and bottom of the hive
lPvTop = HT_SaturationVapourPressure(lSensors(TOP_T).data, 'kelvin') .* lSensors(TOP_RH).data;
lPvBot = HT_SaturationVapourPressure(lSensors(BOT_T).data, 'kelvin') .* lSensors(BOT_RH).data;
lPvExt = HT_SaturationVapourPressure(lSensors(EXT_T).data, 'kelvin') .* lSensors(EXT_RH).data;

lXTick = 0:12:lTimeVec(~isnan(lTimeVec))(end);
lXTickLabel = arrayfun(@(v) sprintf('%.0fH\n%.0f', mod(v, 24), floor(v/24)), lXTick, 'UniformOutput', false);

figure(1);
clf;
subplot(3,1,1);
[ax h1 h2] = plotyy(lTimeVec, lSensors(TOP_RH).data, lTimeVec, lSensors(TOP_T).data-273.15); hold on;
[ax h3 h4] = plotyy(lTimeVec, lSensors(BOT_RH).data, lTimeVec, lSensors(BOT_T).data-273.15);
h5 = plot(ax(2), lTimeVec, lSensors(SOUTH_EXT_T).data-273.15, 'linewidth', 2, 'color', lSensors(SOUTH_EXT_T).color);
h6 = plot(ax(2), lTimeVec, lSensors(EXT_T).data-273.15, 'linewidth', 2, 'color', lSensors(EXT_T).color);
h7 = plot(ax(2), lTimeVec, lSensors(EXT_T2).data-273.15, 'linewidth', 2, 'color', lSensors(EXT_T2).color);
h7 = plot(ax(1), lTimeVec, lSensors(EXT_RH).data, 'linewidth', 2, 'color', lSensors(EXT_RH).color);
lhList = [h1 h2 h3 h4];
set(lhList, 'linewidth', 2);
arrayfun(@(i) set(lhList(i), 'color', lSensors(1+i).color), 1:4);
set(ax, 'xtick', lXTick, 'xticklabel', lXTickLabel, 'xlim', [0 ceil(lXTick(end)/12)*12]);
set(ax(1), 'ycolor', lSensors(4).color);
set(ax(2), 'ycolor', lSensors(5).color);
ylabel(ax(1), 'Humidity [-]', 'fontsize', 16, 'fontweight', 'bold');
ylabel(ax(2), 'Temperature [degC]', 'fontsize', 16, 'fontweight', 'bold');
title('Measurements inside and outside the test hive', 'fontsize', 16, 'fontweight', 'bold');
grid on;

subplot(3,1,2);
plot(lTimeVec, lPvTop, 'color', lSensors(TOP_RH).color, 'linewidth', 2); hold on;
plot(lTimeVec, lPvBot, 'color', lSensors(BOT_RH).color, 'linewidth', 2);
plot(lTimeVec, lPvExt, 'color', lSensors(EXT_RH).color, 'linewidth', 2);


set(gca, 'xtick', lXTick, 'xticklabel', lXTickLabel, 'xlim', [0 ceil(lXTick(end)/12)*12]);
ylabel('Vapour pressure [Pa]' , 'fontsize', 16, 'fontweight', 'bold');
grid on;

subplot(3,1,3);
plot(lTimeVec, lPvTop ./ (8.314/18.01*lSensors(TOP_T).data), 'color', lSensors(TOP_RH).color, 'linewidth', 2); hold on;
plot(lTimeVec, lPvBot ./ (8.314/18.01*lSensors(BOT_T).data), 'color', lSensors(BOT_RH).color, 'linewidth', 2);
plot(lTimeVec, lPvExt ./ (8.314/18.01*lSensors(EXT_T).data), 'color', lSensors(EXT_RH).color, 'linewidth', 2);
xlabel('Time [day]', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Density [kg/m3]' , 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'xtick', lXTick, 'xticklabel', lXTickLabel, 'xlim', [0 ceil(lXTick(end)/12)*12]);
grid on;

% ==============================================================================
% Computation of water transfer through hive walls

% Partial pressure vapour outside the hive
lPvSouthOut = HT_SaturationVapourPressure(lSensors(EXT_T).data, 'kelvin') .* lSensors(EXT_RH).data;
% Relative humidity of the southern outside face at equilibrium
lHRSouthOut = lPvSouthOut ./ HT_SaturationVapourPressure(lSensors(SOUTH_EXT_T).data, 'kelvin');

% Temperature coefficient computed with standard wood equilibrium diagram between HR,T=60%,13°C et 60%,40°C
lTempCoef = -0.00337;
lSorptionf = @(RH, T) 0.244*exp(0.69*log(RH).*exp(1.69*RH)).*(1 + lTempCoef*(T-293.15) );

figure(2);
clf;
% Temperature and humidity inside the hive
subplot(2,2,1);
[ax h1 h2] = plotyy(lTimeVec, lSensors(EXT_RH).data, lTimeVec, lSensors(EXT_T).data-273.15);
set(h1, 'color', lSensors(EXT_RH).color, 'linewidth', 2, 'linestyle', ':', 'displayname', lSensors(EXT_RH).name);
set(h2, 'color', lSensors(EXT_T).color, 'linewidth', 2, 'displayname', lSensors(EXT_T).name);
set(ax, 'fontsize', 16, 'fontweight', 'bold');
set(ax(1), 'ylim', [0 1.1], 'ycolor', 'black');
set(ax(2), 'ylim', [-5 25], 'ycolor', 'black');
##set(get(ax(2), 'ylabel'), 'string'
ylabel(ax(2), 'Temperature [degC]', 'fontsize', 16, 'fontweight', 'bold');
ylabel(ax(1), 'Relative humidity [-]', 'fontsize', 16, 'fontweight', 'bold');
title('Hygrothermal conditions inside the hive', 'fontsize', 16, 'fontweight', 'bold');
set(ax, 'xlim', [10 150]);
grid on;

subplot(2,2,2);
% Temperature and humidity outside the hive
% Relative humidity corrected with the outer face temperature, ie relative humidity
% when the air is coming close to the surface. The surface temperature is measured using
% an infrared camera with special infrared window.
[ax h1 h2] = plotyy(lTimeVec, lSensors(EXT_RH).data, lTimeVec, lSensors(EXT_T).data-273.15); hold on;
[ax h3 h4] = plotyy(lTimeVec, lHRSouthOut, lTimeVec, lSensors(SOUTH_EXT_T).data-273.15);
h5 = plot(ax(2), lTimeVec, lSensors(SOUTH_EXT_PTAT).data-273.15);

set(h1, 'color', lSensors(EXT_RH).color, 'linewidth', 2, 'linestyle', ':', 'displayname', lSensors(EXT_RH).name);
set(h2, 'color', lSensors(EXT_T).color, 'linewidth', 2, 'displayname', lSensors(EXT_T).name);
set(h3, 'color', 'red', 'linewidth', 2, 'linestyle', ':', 'displayname', 'RH surf');
set(h4, 'color', lSensors(SOUTH_EXT_T).color, 'linewidth', 2, 'displayname', lSensors(SOUTH_EXT_T).name);
set(ax, 'fontsize', 16, 'fontweight', 'bold');
ylabel(ax(2), 'Temperature [degC]', 'fontsize', 16, 'fontweight', 'bold');
ylabel(ax(1), 'Relative humidity [-]', 'fontsize', 16, 'fontweight', 'bold');
set(ax, 'xlim', [10 150]);
set(ax(1), 'ylim', [0 1.4], 'ycolor', 'black');
set(ax(2), 'ylim', [-10 60], 'ycolor', 'black');
##title('Hygrothermal conditions outside the hive', 'fontsize', 16, 'fontweight', 'bold');

grid on;
hl = legend();
set(hl, 'numcolumns', 2, ...
        'location', 'north', ...
        'fontsize', 12);
lPos = get(hl, 'position');
set(hl, 'position', lPos + [0 0.08 0 0]);

% Equilibrium Water content
subplot(2,2,3);
lWeqIn  = lSorptionf(lSensors(EXT_RH).data, lSensors(EXT_T).data);
lWeqOut = lSorptionf(min(lHRSouthOut, 1), lSensors(SOUTH_EXT_T).data);

lValid = ~(isna(lWeqIn) | isnan(lWeqIn) | isna(lTimeVec) | isnan(lTimeVec));
lWeqInAverage = trapz(lTimeVec(lValid), lWeqIn(lValid))/(lTimeVec(lValid)(end)-lTimeVec(lValid)(1));
lWeqOutAverage = trapz(lTimeVec(lValid), lWeqOut(lValid))/(lTimeVec(lValid)(end)-lTimeVec(lValid)(1));

plot(lTimeVec, lWeqIn, 'linewidth', 2, 'color', 'black', 'displayname', 'input'); hold on;
plot(lTimeVec, lWeqOut, 'linewidth', 2, 'color', 'red', 'displayname', 'output');
set(gca, 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'xlim', [10 150]);
grid on;
xlabel('Time [hour]', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Water content [-]', 'fontsize', 16, 'fontweight', 'bold');
title(sprintf('Average: WeqIn %.2f ; WeqOut %.2f', lWeqInAverage, lWeqOutAverage), 'fontsize', 16, 'fontweight', 'bold');


tVec = 3600*lTimeVec(lValid);
lWeqInValid = lWeqIn(lValid);
lWeqInFunc = @(t, i, W) interp1(tVec, lWeqInValid, t);

lWeqOutFunc = @(t, i, W) interp1(tVec, lWeqOut(lValid), t);

tVecSimulation = tVec(1)+(0:0.002:1)'*(tVec(end)-tVec(1));

lParamsNLCase = struct(...
                 'length', 0.021, ...                       ... [m] length of the wall
                 'diffusion', @(t,i,W) 2.1E-11*exp(5.77*W), ... [function_handle D = @(W)] diffusion coefficient
                 'convection', [1E-5, 1E-5],                ... [1x2] or [nt x 2] convection coefficient for each sides of the wall
                 'weq', {{lWeqInFunc, lWeqOutFunc}},                  ... [1x2] or [nt x 2] equilibrium water content on each sides
                 'winit', 0.2,                              ... [1x1] or [nx1] initial water content
                 'n', 21, ...                               ... number of nodes
                 'gridType', 'ff', ...
                 't', tVecSimulation);

lOptions = struct('verbose', true);
[Wmat Winfos] = HT_Standalone_Solve_Hygro1D(lParamsNLCase, lOptions);


figure(3);
clf;
nsub = 1:4:numel(tVecSimulation);
surf(tVecSimulation(nsub)/3600/24, Winfos.nodePosition, Wmat(:,nsub));
set(gca, 'fontsize', 16, 'fontweight', 'bold');
zlabel('Water content', 'fontsize', 16, 'fontweight', 'bold');
xlabel('Time [day]', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Normalized position [-]', 'fontsize', 16, 'fontweight', 'bold');
title('Evolution of water content inside hive walls vs time', 'fontsize', 16, 'fontweight', 'bold');

set(gca, 'cameraposition', [4.305830E+01   7.140232E+00   1.0745623E+00]);
set(gca, 'cameratarget', [ 6.000000000000000   0.500000000000000   0.150000000000000]);

