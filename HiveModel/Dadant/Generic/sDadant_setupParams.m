%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022: Montpellier University / CoActions-AltRD-Emmanuel Ruffio
%  Author: emmanuel.ruffio@gmail.com
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
% See Model_Dadant_Schemas.pdf
% ========================================================================

% ========================================================================
%                              PARAMETERS
% ========================================================================

% Definition of materials
lMatWood        = HT_Material_New('wood');
lMatPolystyrene = HT_Material_New('polystyrene');
lMatAir         = HT_Material_New('air');
lMatSteal       = HT_Material_New('steal');

% Global parameters
lHiveParams = struct(...
  'azimuth', 0,                         ... % [radian] Orientation of the hive with respect to geographic north (azimuth=pi=180°=south) (pi/2=east)
  'tilt', 0,                            ... % [radian] Hive tilt
  'dimension', [0.5; 0.5],              ... % [dim2x1] [m] Hive dimensions along X,Y direction
  'bodyHeight', [0.3; 0.2],             ... % [m] Height of hive body part 1 and 2
  'wallThickness', 0.02,                ... % [m] Thickness of hive walls
  'overframeThickness', 5E-3,           ... % [m] Thickness of overframe layer
  'underroofThickness', 1E-2,           ... % [m] Thickness of the layer between overframe and roof
  'roofThickness', 0.6E-3,              ... % [m] Thickness of hive roof
  'roofIntExtra', [50E-3; 50E-3],       ... % [dim2x1] [m] Extra space inserted between roof and sidewalls (along X,Y dimension). Thickness of air layer.
  'roofHeight', [0.10],                 ... % [m] Height of roof sides
  'longitude', [43 36 39],              ... % [deg minute second] Longitude
  'latitude', [3 52 38],                ... % [deg minute second] Latitude
  'date', [2023 01 10],                 ... % [year month day] Date
  'time', [0 0 0],                      ... % [hour minute second] Time
  'globalAxis', eye(3),                 ... % Identity matrix denotes the 3 main direction of X,Y,Z in columns
  'globalPosition', zeros(3,1)          ... % [m x m x m] Position of the hive
  );

lHiveMat = struct(...                   ... % Structure containing hive materials
  'sidewalls_int',  lMatWood,           ... % Material used to specify thermal emissivity inside the hive (it is not used for the side walls material)
  'sidewalls',      lMatWood,           ...
  'overframe',      lMatWood,           ... % Used for bulk material parameters of high side (toward Z+) emissivity
  'overframe_int',  lMatWood,           ... % Used only to define the lower side emissivity
  'underroof',      lMatAir,            ... % Layer between overframe and the roof
  'roof_int',       lMatSteal,          ... % Distinguish internal and external is useful to define different emissivities
  'roof_ext',       lMatSteal,          ...
  'roofsides_int',  lMatSteal,          ...
  'roofsides_ext',  lMatSteal,          ...
  'roofsideAir',    lMatAir,            ...
  'externalAir',    lMatAir,            ...
  'bodyAir',        lMatAir
  );

lGroundParams = struct(...
  'material', HT_Material_New('wet ground')    ... % Material of the ground (not used in this model yet)
  );

lConvection = struct(...
  'walls1_int', { {5, 5, 5, 5, 5, NaN} },                   ... % [W/m²/K] Convection coefficient for side walls 1 and bottom (internal)
  'walls1_ext', { {10, 10, 10, 10, 10, NaN} },              ... % [W/m²/K] Convection coefficient for side walls 1 and bottom (external xm,xp,ym,yp)
  'walls2_int', { {5, 5, 5, 5, NaN, 5} },                   ... % [W/m²/K] Convection coefficient for side walls 2 and top (internal)
  'walls2_ext', { {1, 1, 1, 1, NaN, NaN} },                 ... % [W/m²/K] Convection coefficient for side walls 2 (external)
  'roof_ext',   { {10, 10, 10, 10, NaN, 10} }   ... % [W/m²/K] Convection coefficient for roof (sides and top)
  );

lContact = struct(...                               % Structure containing information about contact resistances
  'overframe_roof', 0.5E-3/lMatAir.lambda       ... % [Km²/W] Contact resistance between hive roof and wood below
  );

lCommand = struct(...
  'type', 'simulation', ...
  'scriptFile', 'Dadant_setupCommandsDefault');

lRadiation = struct(...
  'enable_hiveBody', false,                     ... % [logical] enable/disable radiation heat flow inside hive body
  'enable_overframe', false                     ... % [logical] enable/disable radiation heat flow inside overframe air layer
  );

lMesh = struct(...
  'hivebody', [5 5],                            ... % Number of nodes in hive part1 and 2
  'sidewalls', 5,                               ... % number of nodes inside each hive side walls 1 and 2, bottom
  'overframe', [3 3 4],                         ... % number of nodes inside overframe layer in X,Y,Z direction
  'underroof', 4,                               ... % Number of nodes in z-axis direction of the underroof layer (other directions are taken from .overframe)
  'roof', 1,                                    ... % Number of nodes of roof (in depth(=z) axis)
  'airside', 2                                  ... % Number of nodes of air side layer (in depth axis)
  );

lComputation = struct(...
  'nt',         1000,                           ... % Number of time steps
  'startTime',  0,                              ... % [second] Start time
  'timeStep', 24*3600/1000,                     ... % [second] Time step duration
  'initTemperature', 10 ...25                         ... % [degres] Initial temperature
  );

lOptions = struct('verbose', true,      ... % Display additional information
                  'verboselevel', 4,    ... % Display filter (the higher the more message)
                  'checkGeometry', true ... % Perform additional tests to check geometry coherency
                  );

l3DPlotOptions = struct(...
  'explode', 0, ...
  'explodeCenter', [0.5*lHiveParams.dimension; 0.5*sum(lHiveParams.bodyHeight)], ...
  'displaygrid', true, ...
  'displayedge', false, ...
  'facecolortype', 'flat'...
  );

lPlotOptions = struct(...
  'onlyOneFigure', false, ...
  'plotGeometry', true, ...                % Plot the geometry in figure 1
  'timeStepPolicy', 'none', ...             % 'none' : keep the time step equal to computation time step
                                            % 'fixed' : change time step according to the field "fixedTimeStep"
  'fixedTimeStep', 100, ...
  'axisFontSize', 16, ...
  'axisFontWeight', 'bold', ...
  'labelFontSize', 16, ...
  'labelFontWeight', 'bold', ...
  'titleFontSize', 16, ...
  'titleFontWeight', 'bold' ...
  );

lPlotColors = struct(...
  'faces', [ 0.00000   0.44700   0.74100; ... XM % Display colors used to plot temperature for each face
             0.85000   0.32500   0.09800; ... XP % It allows faces to have always the same color
             0.92900   0.69400   0.12500; ... YM
             0.49400   0.18400   0.55600; ... YP
             0.46600   0.67400   0.18800; ... ZM
             0.30100   0.74500   0.93300; ... ZP
             0.63500   0.07800   0.184000],... Roof
  'hbody', struct('bottom', [0  1  0], 'top', [0.5 1 0.5]), ...
  'roof', [0.63500   0.07800   0.184000]);

lSensitivityParams = struct(...
  'multEpsilon', 0.02, ...
  'backupFile', 'HiveModel_Dadant_SensitivityBackup.mat', ...
  'backupReset', true);

lSensitivityAnalysis = struct(...
  'name', 'wall lambda',...
  'fields', {{'hiveMat', 'sidewalls'}}, ...
  'setterNext', @(params) setfield(params, {1}, 'hiveMat', {1}, 'sidewalls', HT_Material_SetLambda(params.hiveMat.sidewalls, 'multiply', 1+lSensitivityParams.multEpsilon)), ...
  'setterPrev', @(params) setfield(params, {1}, 'hiveMat', {1}, 'sidewalls', HT_Material_SetLambda(params.hiveMat.sidewalls, 'multiply', 1-lSensitivityParams.multEpsilon)), ...
  'getUserData', @(params) HT_Material_GetLambda(params.hiveMat.sidewalls, 1), ...
  'compute', @(userData, Xp, Xn) (Xn-Xp) / (2*userData*lSensitivityParams.multEpsilon));

lSensitivityAnalysis = [lSensitivityAnalysis;
  struct(...
  'name', 'overframe lambda',...
  'fields', {{'hiveMat', 'overframe'}}, ...
  'setterNext', @(params) setfield(params, {1}, 'hiveMat', {1}, 'overframe', HT_Material_SetLambda(params.hiveMat.overframe, 'multiply', 1+lSensitivityParams.multEpsilon)), ...
  'setterPrev', @(params) setfield(params, {1}, 'hiveMat', {1}, 'overframe', HT_Material_SetLambda(params.hiveMat.overframe, 'multiply', 1-lSensitivityParams.multEpsilon)), ...
  'getUserData', @(params) HT_Material_GetLambda(params.hiveMat.sidewalls, 1), ...
  'compute', @(userData, Xp, Xn) (Xn-Xp) / (2*userData*lSensitivityParams.multEpsilon));
  ];

lParams = struct();
lParams.hiveParams = lHiveParams;
lParams.hiveMat = lHiveMat;
lParams.groundParams = lGroundParams;
lParams.convection = lConvection;
lParams.contact = lContact;
lParams.command = lCommand;
lParams.radiation = lRadiation;
lParams.mesh = lMesh;
lParams.computation = lComputation;
lParams.plotOptions = lPlotOptions;
lParams.plotColors = lPlotColors;
lParams.plotOptions3d = l3DPlotOptions;
lParams.options = lOptions;





