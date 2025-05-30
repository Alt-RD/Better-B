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
% Example: using library HiveTemp.
% This script shows the heat radiation inside a empty box. The bottom face
% is then heated with a constant temperature and the other face temperatures
% rise slowly due to radiation heat transfer. Walls are made of wood.
%
% See: HiveTemp Manual.pdf for more details.
% ========================================================================

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../';
addpath(HT_VAR_LIB_PATH);

% Init the library
HT_Init();

% Definition of materials
lMatWood = HT_Material_New('wood', 'eps',    [HT_WAVELENGTH_VISIBLE  0.5; ...
                                              HT_WAVELENGTH_IR_MED   1]);
lMatAir  = HT_Material_New('air');

% Defines the global axis (orientation of the hive) and position
lGlobalAxis = eye(3);              % Identity matrix denotes the 3 main direction of X,Y,Z in columns
lGlobalPosition = [0.0; 0.0; 0.0]; % [m x m x m] Position of the hive

% Global parameters
lOptions = struct('verbose', true);  % Display additional information (level 4)
lCheckGeometry = true;            % Perform additional tests to check geometry coherency
                                  % (in HT_Model_ConnectFaces function)

lBodySize = [0.5 0.5 0.02]';       % [m x m x m]: Hive body size along X,Y,Z axis direction
lWallThickness = 0.01;            % [m]: Thickness of the box walls

lNodeCount = 2;                   % Number of wall nodes along the length-axis.

lInitTemperature = 25;            % [°C or K] : Initial temperature
lCmdTemperature = 100;

tSpec = {0 30*3600/1000 100};    % [s] : Time vector specification {t0, tf, dt}

l3DPlotOptions = struct(...
  'explode', 0, ...
  'explodeCenter', 0.5*lBodySize, ...
  'displaygrid', true, ...
  'displayedge', true, ...
  'facecolortype', 'flat',...
  'edgecolor', [0 0 0]...
  );

lDisplayOptions = struct( 'explode', 1.2, ...
                          'explodeCenter', 0.5*lBodySize, ...
                          'displaygrid', true, ...
                          'displayedge', false, ...
                          'facecolortype', 'flat');

figure(1);
clf;

% ========================================================================
%                               MODEL
% ========================================================================
% Create model main object with name "main"
lModel = HT_Model_Init('', 'main');

% Create an empty box.
% The model is discarded. We are only interested in 6 outside faces
[~, lMod_Box_faces] = HT_Model_Conduction1D(                       ... % Build a two layers 1D heat conduction problem
                          "mEmpty",                                     ... % Model name of the hive body part 1 and 2
                          'material', lMatAir,                          ... % Air material for the hive body parts
                          'dim', lBodySize([X Y]),                      ... % Dimensions along X and Y direction
                          'axis', lGlobalAxis,                          ... % Local axis is equal to global axis
                          'length', lBodySize(Z),... % Length of the two layers (part1 and part2)
                          'n', 1,                                       ... % Number of nodes in each layer
                          'mergeFaces', false,...
                          'gridType', 'ff', ...
                          lOptions);

% This cell array will store the internal faces of the box
% This faces differ from lMod_Box_faces since they have a mesh defined
% but the geometry is the same.
lIntFaces = cell(6,1);

for i=1:6
  [lMod_SideWall, lMod_SideWall_faces] = HT_Model_Conduction3D(...
    sprintf('mSW_%s', FaceDIRstrL{i}), ...
    struct( 'base',     lMod_Box_faces(i),          ... % Retrieve the face of the empty box to build a volume based on it
            'length',   lWallThickness,             ... % Tickness of the walls
            'material', lMatWood,                   ...
            'n',        [lNodeCount lNodeCount 1],  ... % NA is used since mesh is extracted from face <base>
            'zGridType', 'ff',                       ... % 'f' defines the size of boundary nodes ('h' for halved, or 'f' for full)
            'xGridType', 'ff', ...
            'yGridType', 'ff'
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_SideWall, lOptions);

  lIntFaces{i} = lMod_SideWall_faces(FaceZM);

  HT_Plot_Face(lIntFaces{i}, l3DPlotOptions); % Display the geometry
endfor

% Create the radiation model and retrieve the fluxMatrix [-]
% For a given vector of temperature T, the amount of incoming radiation
% for each nodes is given by Phi[W/m²] = lFluxMatrix * sigma * T^4
% HT_SolveModel does the same but multiply by node areas to get the absolute
% heat flux [W]
[lModel_Radiation lFluxMatrix] = HT_Model_RadiationInBox('mRadBox', ...
                                            struct( 'face', { lIntFaces }, ...
                                                    'wavelength', HT_WAVELENGTH_IR_MED));
lModel = HT_Model_Merge(lModel, lModel_Radiation, lOptions);

% ========================================================================
%                               COMMAND
% ========================================================================
lCmd = HT_Cmd_Init(
          'name', 'cmd_bottom', ...
          'type', 'temperature', ...
          'nodes', lIntFaces{FaceXM}, ... Nodes will be extracted from this face
          'data', lCmdTemperature);

% ========================================================================
%                             RESOLUTION
% ========================================================================
[T, ~, Tnodes, ~, t] = HT_SolveModel(...
         lModel,                                            ... Model to be solved
         lCmd,                                              ... List of commands
         lInitTemperature,                                  ... Initial temperature
         tSpec,                                             ... Time vector
         struct('all', true,                                ... All temperature must be returned. T will be (dim Nnodes x Ntimes)
                'algorithm', 'linear',                      ... Algorithm used
                'replaceT0', true,                          ... Overwrite all initial temperature with specified temperature <lInitTemperature>
                'verbose', lOptions.verbose,                ...
                'unit', 'degres',                           ... % Warns the function that all temperature are expressed in degC.
                'info', 'default'));

% ========================================================================
%                             RESULTS
% ========================================================================

T_XPface = HT_Result_GetFaceTemperature(lIntFaces{FaceXP}, T, Tnodes, 'timevector', t, 'operation', 'average');
T_YMface = HT_Result_GetFaceTemperature(lIntFaces{FaceYM}, T, Tnodes, 'timevector', t, 'operation', 'average');
T_XPnodes = HT_Result_GetFaceTemperature(lIntFaces{FaceXP}, T, Tnodes, 'timevector', t, 'operation', 'array');

figure(1);
xlabel('X', 'fontsize', 20, 'fontweight', 'bold');
ylabel('Y', 'fontsize', 20, 'fontweight', 'bold');
zlabel('Z', 'fontsize', 20, 'fontweight', 'bold');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
grid on;
camtarget([0.01 0.25 0.25]);
campos([-0.1042  -2.1648 2.4343]);

figure(2);
clf;
plot(t, T_XPface, '-', 'displayname', 'Average XP face'); hold on;
plot(t, T_YMface, '-', 'displayname', 'Average YM face'); hold on;
for i=1:rows(T_XPnodes)
  plot(t, T_XPnodes(i,:), 'displayname', sprintf('XPface.node%d', i));
endfor
grid on;
hl = legend();
set(hl, 'fontsize', 14);
title('Temperature evolution', 'fontsize', 20, 'fontweight', 'bold');


