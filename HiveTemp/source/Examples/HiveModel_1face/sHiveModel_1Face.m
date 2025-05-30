%  This file is part of project HiveTemp.
%
%  Copyright (c) CoActions-AltRD-Emmanuel Ruffio
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
% Programme d'exemple de l'utilisation de la librairie HiveTemp.
% Ce script construit un modèle très simplifié d'une ruche.
% Un volume d'air est séparé de l'environnement par un mur en bois d'épaisseur
% 2cm. L'environnement est à 10°C. La température initiale de la ruche est de 10°C.
%
% ========================================================================

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

lGlobalAxis = eye(3);              % Identity matrix denotes the 3 main direction of X,Y,Z in columns
lGlobalPosition = [0.0; 0.0; 0.0]; % [m x m x m] Position of the hive

% Global parameters
lOptions = struct('verbose',      true, ... % Display additional information (level 4)
                  'verboseLevel', 1);
lCheckGeometry = true;                      % Perform additional tests to check geometry coherency
                                            % (in HT_Model_ConnectFaces function)

lBodySize = [0.5 0.5 0.5]';   % [m x m x m]: Hive body size along X,Y,Z axis direction
lInternalConvectionCoef = 5;  % [W/m²/K] : COnvection coefficient between internal air (hive body)
                              % and surroundings hive walls
lExternalConvectionCoef = 10;  % [W/m²/K] : COnvection coefficient between internal air (hive body)
                              % and surroundings hive walls
lWallLength = 0.02;           % [cm] : Hive wall thickness
lNodeCount = 10;              % Number of wall nodes along the length-axis.

lInitTemperature = 20;        % [°C or K] : Initial temperature


% View parameters
lViewOptions = struct('explode', 0.0);

% Definition of materials
% Custom materials
lMatwood = HT_Material_New( 'name', 'bois', ...
                            'rho', 600, ...
                            'cp', 1600, ...
                            'lambda', 0.16, ...
                            'eps', 0.9);

lMatpolystyrene = HT_Material_New( 'name', 'polystyrene', ...
                            'rho', 20, ...
                            'cp', 1450, ...
                            'lambda', 0.035, ...
                            'eps', 0.9);

% Or builtin material
lMatair = HT_Material_New( 'air' );

% Create model main object with name "main"
lModel = HT_Model_Init('', 'main');

% Create the hive body with dimension 0.5x0.5x0.5m
[lModel0D lModel0D_faces] = HT_Model_Conduction0D("model_hiveBody", ...
                  struct('dim', lBodySize, ...
                         'material', lMatair, ...
                         'axis', lGlobalAxis, ...
                         'globalPosition', lGlobalPosition), ...
                  setfield(lOptions, 'nodeName', 'nbody')); % Specify the name of the node (instead of default name)
% Merge the 0D Model to the main model
lModel = HT_Model_Merge(lModel, lModel0D, lOptions);

% Create a 1-dimensionnal single layer wall model in the XP (X+) direction
[lModelWall1D lModelWall1D_faces lModelWall1D_nodes] = HT_Model_Conduction1D("model_hiveWall", ...
                  'length', lWallLength,                  ... % Thickness of the wall
                  'material', lMatwood,                   ... % Material
                  'base', lModel0D_faces(FaceXP),         ... % Dimensions of the wall are equal to the body size along Y and Z
                  'n', lNodeCount,                        ... % Spatial discretization. Wall is divided into "n" nodes
                  'boundaryNames', {{'nwallxp_int', 'nwallxp_out'}}, ... % Double {{ }} are necessary since there is a underlying structure initialization. {} would mean structure array.  That's not what we want
                  'gridType', 'hh', ...
                  lOptions);
% Merge the 1D Model to the main model
lModel = HT_Model_Merge(lModel, lModelWall1D, lOptions);
clear lModelWall1D;

% Create a thermal connection between hive body and surrounding wall
% Connection with wall XP
lModel_BodyWallXm = HT_Model_ConnectFaces('model_body-wallxm', ...
                        lModel0D_faces(FaceXP),  ...
                        lModelWall1D_faces(FaceXM), ...
                        struct('g', lInternalConvectionCoef),
                        setfield(lOptions, "checkGeometry", lCheckGeometry)); % Add a field to the structure lOptions
% Merge the model of the thermal connection to the main model
lModel = HT_Model_Merge(lModel, lModel_BodyWallXm, lOptions);
clear lModel_BodyWallXm;

lAirExtNode = HT_Node_Init('AirExt');
lNodeModel = HT_Model_Connect('model_AirExtConnection', ...
                              lAirExtNode, ...
                              lModelWall1D_faces(FaceXP), ... % Retrieve the external face
                              'g', lExternalConvectionCoef, ... Uniform conductance between external face and outside air
                              lOptions);
lModel = HT_Model_Merge(lModel, lNodeModel, lOptions);
clear lAirExtNode;

% ======================= Command definition ==========================
% Retrieve the area of external hive
lWall1D_extFaceArea = HT_Face_GetAbsoluteData(lModelWall1D_faces(FaceXP), 'nodesArea');

lCmd = HT_Cmd_Init(
        'name', 'airext', ...
        'type', 'temperature', ...
        'nodes', 'AirExt', ...
        'data', 10, ... External temperature of 10°C

        'name', 'sun', ...
        'type', 'flux', ...
        'nodes', {'nwallxp_out'}, ...
        'area', lWall1D_extFaceArea, ...
        'data', 1000); ... Solar irradiance of 1000W/m²

% ======================= Model resolution ==========================
t = (0:0.001:1)*4*3600; % [s] : Time vector

[T, U, Tnodes] = HT_SolveModel(...
         lModel, ...  Model to be solved
         lCmd, ...    List of commands
         lInitTemperature, ... Initial temperature
         t, ... Time vector
         struct('all', true,...
                'algorithm', 'chol',...
                'replaceT0', true, ...
                'verbose', lOptions.verbose));

% ======================= Display geometry ==========================
figure(1);
clf;

% Display the 6 faces of the 0D Model
HT_Plot_Face(lModel0D_faces);
% Display the 6 faces (+ eventually the spatial discretization) of the XP wall
HT_Plot_Face(lModelWall1D_faces);
grid on;
xlabel('X', 'fontsize', 14, 'fontweight', 'bold');
ylabel('Y', 'fontsize', 14, 'fontweight', 'bold');
zlabel('Z', 'fontsize', 14, 'fontweight', 'bold');

% ======================= Display temperature evolution ========================
figure(2);
clf;
% Extract the temperature evolution of node called "nbody"
tplot = 0:100:t(end);
Thivebody = HT_Result_GetNodeTemperature('nbody', T, Tnodes, ...
                                         'time', tplot, 'timevector', t);
plot(tplot/3600, Thivebody, '-+', 'linewidth', 2.0); hold on;
set(gca, 'fontsize', 14);
xlabel('Temps [hour]');
ylabel('Temperature [degC]');
title('Temperature evolution of hive body temperature');
grid on;

% ================ Display temperature profils inside Xp wall ==================
figure(3);
clf;
% Extract temperatures for each nodes of XP wall. The list of nodes was saved
% in <lModelWall1D_nodes>. The time is interpolated and is uniformly-spaced.
tprofil = logspace(0, log10(t(end)), 10); %0:600:t(end); % [s]
nNodes = numel(lModelWall1D_nodes);
Twallprofil = HT_Result_GetNodeTemperature(lModelWall1D_nodes, T, Tnodes, ...
                                          'time', tprofil, ... Time of different profils
                                          'timevector', t);
for k=1:size(Twallprofil,2)
  plot(Twallprofil(:,k), 'linewidth', 2.0); hold on;
  text(nNodes/2, Twallprofil(idivide(nNodes,int32(2)),k), ...
        sprintf('t=%.1fh', tprofil(k)/3600),...
        'fontsize', 14, ...
        'fontweight', 'bold');
endfor

set(gca, 'fontsize', 14);
xlabel('Position [node]');
ylabel('Temperature [degC]');
title('Temperature profils inside hive wall');
grid on;

