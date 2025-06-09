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
%
% See Model_TruncHive_Schemas.pdf
% ========================================================================

% ========================================================================
%                              PARAMETERS
% ========================================================================

% Definition of materials
lMatWood        = HT_Material_New('wood');
lMatAir         = HT_Material_New('air');
lMatGround      = HT_Material_New('dry ground'); % https://energieplus-lesite.be/donnees/enveloppe44/caracteristiques-thermiques-des-sols/
lMatInsulation  = HT_Material_New('polystyrene', 'lambda', 0.05);

% Global parameters
lHiveParams = struct(...
  'azimuth', 0,                         ... % [radian] Orientation of the hive with respect to geographic north (azimuth=pi=180°=south) (pi/2=east)
  'extRadius', 0.4,                     ... % [m] External Hive radius of body part computed between two parallel sides
  'bodyHeight', 0.8,                    ... % [m] Height of hive body part 1 and stone
  'ovfHeight', 0.15,                     ... % [m] Minimal height of the overframe part
  'roofTilt', 20/360*2*pi,              ... % [radian] Title angle of the roof
  'wall', struct(...
    'intThickness', 0.015, ...          ... % [m] Thickness of internal layer of hive walls
    'extThickness', 0.015, ...          ... % [m] Thickness of external layer of hive walls
    'insulationThickness', 0.015,       ... % [m] Thickness of insulation layer of hive walls
    'baseThickness', 0.025             ... % [m] Thickness of the hive base
    ),...
  'overframeThickness', 0.025,          ... % [m] Thickness of roof wood layer
  'roofThickness', 0.025,               ... % [m] Thickness of hive roof
  'ovfContentVolume', 0.2*0.2*0.02,     ... % [m3] Volume of wood in the overframe air part
  'longitude', [43 36 39],              ... % [deg minute second] Longitude
  'latitude', [3 52 38],                ... % [deg minute second] Latitude
  'date', [2025 01 10],                 ... % [year month day] Date
  'time', [0 0 0],                      ... % [hour minute second] Time
  'globalAxis', eye(3),                 ... % Identity matrix denotes the 3 main direction of X,Y,Z in columns
  'globalPosition', zeros(3,1)          ... % [m x m x m] Position of the hive
  );

% By default, it is assumed that first basis vector is oriented towards south.

lHiveMat = struct(...                   ... % Structure containing hive materials
  'walls',          lMatWood,           ...
  'extWalls',       lMatWood,           ... % Used (for outside emissivity)
  'wallInsulation', lMatInsulation,     ...
  'overframe',      lMatWood,           ... % Used for bulk material parameters of high side (toward Z+) emissivity
  'roof',           lMatWood,           ...
  'internalAir',    lMatAir,            ...
  'externalAir',    lMatAir             ...
  );

lGroundParams = struct(...
  'material', HT_Material_New('wet ground')    ... % Material of the ground (not used in this model yet)
  );

lConvection = struct(...
  'body_bot', 5,                  ... % [W/m²/K] Convection coefficient for bottom faces
  'body_top', 5,                  ... % [W/m²/K] Convection coefficient for bottom faces
  'body_walls', 5,                     ... % [W/m²/K] Convection coefficient for internal faces (sides)
  'ext_walls', 10,                    ... % [W/m²/K] Convection coefficient for external faces
  'ext_bot', 10,                    ... % [W/m²/K] Convection coefficient for external faces
  'ovf_top', { 10 },                  ...
  'ovf_bot', { 5 },                  ...
  'ovf_int', { 5 },                  ...
  'roof_top', { 10 }                 ... % [W/m²/K] Convection coefficient for roof (sides and top)
  );

lCommand = struct(...
  'type', 'simulation', ...
  'scriptfile', 'sEcological_setupCommands');

lMesh = struct(...
  'wall', struct(...
    'intLayer', 3, ...                              ... % Number of nodes of internal layer of hive walls
    'extLayer', 3, ...                              ... % Number of nodes of external layer of hive walls
    'insulationLayer', 3,                           ... % Number of nodes of insulation layer of hive walls
    'bottom', 5,                                    ... % Number of nodes in the base
    'ovf', 5,                                       ... % Number of nodes in the overframe
    'roof', 5),                                     ... % Number of nodes in the roof
  'bodyHeight', 5,                                  ... % Number of nodes along the vertical axis for the internal air part
  'ovfHeight', 3,                                   ... % Number of nodes in the overframe part (before tilt roof)
  'innerRadius', 1                                  ... % Number of nodes along the radius in the air part
  );

lComputation = struct(...
  'nt',         500,                           ... % Number of time steps
  'timeStep', 24*3600/1000,                     ... % [second] Time step duration
  'initTemperature', 10                         ... % [degres] Initial temperature
  );

lOptions = struct('verbose', true,      ... % Display additional information
                  'verboselevel', 4,    ... % Display filter (the higher the more message)
                  'checkGeometry', true ... % Perform additional tests to check geometry coherency
                  );

l3DPlotOptions = struct(...
  'explode', 0, ...
  'explodeCenter', [0; 0; 0.5*lHiveParams.bodyHeight], ...
  'displaygrid', true, ...
  'displayedge', false, ...
  'facecolortype', 'flat', ...
  'edgecolor', [0 0 0],...
  'displayNormal', false,...
  'displayAxis', false...
  );

lPlotOptions = struct(...
  'onlyOneFigure', true, ...
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

lUserPlotOptions = struct(...
  'figureIndex', [2 3 4 5], ...
  'currentFigureIndex', 1, ...
  'plotSunSideDirection', true, ...     % Plot the angle between each hive sides and the sun
  'plotSunPosition', false ...
  );

lPlotColors = struct(...
  'hbody', struct('bottom', [0  1  0], 'top', [0.5 1 0.5]), ...
  'roof', [0.63500   0.07800   0.184000]);

% Check parameters
assert(lHiveParams.extRadius >= sum([lHiveParams.wall.intThickness + lHiveParams.wall.insulationThickness + lHiveParams.wall.extThickness]));
