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
lMatPolystyrene = HT_Material_New('polystyrene');
lMatAir         = HT_Material_New('air');
lMatGranit       = HT_Material_New('granit');
lMatSteal       = HT_Material_New('steal');
lMatGround      = HT_Material_New('dry ground'); % https://energieplus-lesite.be/donnees/enveloppe44/caracteristiques-thermiques-des-sols/
% Argile/limon sec  lambda  = 0.4   0.5  1.0    rhoC = 1.5 – 1.6
% Argile/limon saturé d’eau = 0.9   1.7   2.3   rhoC = 1.6 – 3.4
lMatBrick       = HT_Material_New('brick');

% Global parameters
lHiveParams = struct(...
  'azimuth', 0,                         ... % [radian] Orientation of the hive with respect to geographic north (azimuth=pi=180°=south) (pi/2=east)
  'tilt', 0,                            ... % [radian] Hive tilt
  'bodyRadius', 0.4,                    ... % [m] [m] External Hive radius of body part and stone above
  'bodyHeight', 0.8,                    ... % [m] Height of hive body part 1 and stone
  'baseRadius', 0.6,                    ... % [m] Radius of the brick base
  'baseThickness', 0.07,                ... % [m] Thickness of the brick base
  'groundRadius', 1.5,                  ... % [m] Radius of the ground
  'groundThickness', 0.2,               ... % [m] Thickness of the ground layer
  'stoneThickness', 0.05,               ... % [m] Thickness of the stone roof
  'stoneRadius', 0.8,                   ... % [m] Radius of the stone roof
  'wallThickness', 0.05,                ... % [m] Thickness of hive walls
  'overframeThickness', 0.03,           ... % [m] Thickness of roof wood layer
  'underroofThickness', 0.01,           ... % [m] Thickness of the layer between overframe and roof
  'roofThickness', 0.6E-3,              ... % [m] Thickness of hive roof
  'roofShift', [0 0],                   ... % [m] Move the roof stone towards [x,y] direction (x is south). Warning, the thermal model is not changed by this parameter, only the illumination computation
  'longitude', [43 36 39],              ... % [deg minute second] Longitude
  'latitude', [3 52 38],                ... % [deg minute second] Latitude
  'date', [2023 01 10],                 ... % [year month day] Date
  'time', [10 0 0],                      ... % [hour minute second] Time
  'globalAxis', eye(3),                 ... % Identity matrix denotes the 3 main direction of X,Y,Z in columns
  'globalPosition', zeros(3,1)          ... % [m x m x m] Position of the hive
  );

% By default, it is assumed that first basis vector is oriented towards south.

lHiveMat = struct(...                   ... % Structure containing hive materials
  'walls_int',      lMatWood,           ... % Material used to specify thermal emissivity inside the hive (it is not used for the side walls material)
  'walls',          HT_Material_SetEmissivity(lMatWood, 'model', 'constant', 'parameters', 4),          ...
  'overframe',      lMatWood,           ... % Used for bulk material parameters of high side (toward Z+) emissivity
  'overframe_int',  lMatWood,           ... % Used only to define the lower side emissivity
  'underroof',      lMatAir,            ... % Layer between overframe and the roof
  'roof',           lMatGranit,          ...
  'internalAir',    lMatAir,            ...
  'externalAir',    lMatAir,            ...
  'ground',         lMatGround,         ...
  'base',           lMatBrick           ... % Used between the ground and the trunc hive
  );

lGroundParams = struct(...
  'material', HT_Material_New('wet ground')    ... % Material of the ground (not used in this model yet)
  );

lConvection = struct(...
  'walls_bottom', 5,                  ... % [W/m²/K] Convection coefficient for bottom faces
  'walls_int', 5,                     ... % [W/m²/K] Convection coefficient for internal faces (sides)
  'walls_ext', 10,                    ... % [W/m²/K] Convection coefficient for external faces
  'ovf_top', { 10 },                  ...
  'ovf_bot', { 10 },                  ...
  'roof_top',   { 10 },               ... % [W/m²/K] Convection coefficient for roof (sides and top)
  'roof_bot',   { 10 },               ...
  'roof_ext',   { 10 }, ...
  'roof_int',   { 10 }, ...
  'base_top', 10,                     ... % [W/m²/K] Convection coefficient between brick base top part and the external air
  'base_side', 10,                    ... % [W/m²/K] Convection coefficient between brick base sides and the external air
  'ground', 10
  );

lCommand = struct(...
  'type', 'simulation', ...
  'scriptfile', 'sHiveModel_Trunk_setupCommands');

lMesh = struct(...
  'bodyRadius', 3,                                  ... % Number of nodes in hive part1 and 2
  'bodyHeight', 6,                                  ... % Number of nodes along the vertical axis
  'bodyAirHeight', 3,                               ... % Number of nodes along the vertical axis for the internal air part
  'innerRadius', 1,                                 ... % Number of nodes along the radius in the air part
  'circular', 8,
  'overframe', 2,                                   ... % number of nodes inside overframe layer in Z direction
  'underroof', 1,                                   ... % Number of nodes in z-axis direction of the underroof layer (other directions are taken from .overframe)
  'roofInner', 1,                                   ... % Number of nodes of roof in the inner cylinder part.
  'roofOuter', 3,                                   ... % Number of nodes of roof in the outer cylinder part (> hiveRadius + wallThickness)
  'roofHeight', 2,                                  ... % Number of nodes of roof along the vertical direction
  'baseOuter', 1,                                   ... % Number of nodes of the base outside the trunc hive
  'baseHeight', 3,                                  ... % Number of nodes of the base along the vertical direction
  'groundHeight', 3,                                ... % Number of nodes in the ground along the vertical direction
  'groundOuter', 3,                                 ... % Number of nodes in the groud outside the brick radius
  'groundOuterRatio', 2                             ... % For >= 1, condense nodes towards the center
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
  'explode', 1, ...
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
  'plotGeometry', false, ...                % Plot the geometry in figure 1
  'timeStepPolicy', 'none', ...             % 'none' : keep the time step equal to computation time step
                                            % 'fixed' : change time step according to the field "fixedTimeStep"
  'fixedTimeStep', 100 ...
  );

lPlotColors = struct(...
  'hbody', struct('bottom', [0  1  0], 'top', [0.5 1 0.5]), ...
  'roof', [0.63500   0.07800   0.184000]);

% Check parameters
assert(lMesh.bodyAirHeight < lMesh.bodyHeight, 'The number of nodes in the air volume must be lower than the number of nodes in the cylinder');
