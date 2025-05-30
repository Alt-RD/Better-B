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
% This script implements a simplified thermal model of a beehive. It is modeled
% as an empty wood box. Thickness of surrounding walls is 2cm. Outside air is
% varying with sinus function. Mean value is 10°C, peak-peak amplitude is 20°C and
% period is 2h.
% Hive initial temperature is 10°C.
%
% Ce script modélise de manière simplifiée le comportement thermique d'une
% ruche. Une ruche est représentée par une boite cubique en bois. Les parois
% latérales ont une épaisseur de 2cm. La température extérieure varie selon une
% sinusoide de valeur moyenne 10°C, d'amplitude crête-crête 20°C et de période
% d'oscillation de 2h. La température initiale de la ruche est de 10°C.
%
% See: HiveTemp Manual.pdf for more details.
% ========================================================================

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

% Definition of materials
lMatwood = HT_Material_New( 'name', 'bois', ...
                            'rho',    600, ...
                            'cp',     1600, ...
                            'lambda', 0.16, ...
                            'eps',    0.9);

lMatpolystyrene = HT_Material_New( 'name', 'polystyrene', ...
                            'rho',    20, ...
                            'cp',     1450, ...
                            'lambda', 0.035, ...
                            'eps',    0.9);

lMatair = HT_Material_New( 'name', 'air', ...
                            'rho',    1.2, ...
                            'cp',     1.009E3, ...
                            'lambda', 0.025);

lMatsteal = HT_Material_New('name', 'steal', ...
                            'rho',    7800, ...
                            'cp',     450, ...
                            'lambda', 50, ...
                            'eps',    0.45);

% Defines the global axis (orientation of the hive) and position
lGlobalAxis = eye(3);              % Identity matrix denotes the 3 main direction of X,Y,Z in columns
lGlobalPosition = [0.0; 0.0; 0.0]; % [m x m x m] Position of the hive

% Global parameters
lOptions = struct('verbose', true);  % Display additional information (level 4)
lCheckGeometry = true;            % Perform additional tests to check geometry coherency
                                  % (in HT_Model_ConnectFaces function)

lHiveAzimuth = 0;                 % Hive orientation with respect to geographic north (azimuth=pi=180°=south) (pi/2=east)
lHiveTilt = 0;                    %
lBodySize = [0.5 0.5 0.5]';       % [m x m x m]: Hive body size along X,Y,Z axis direction

lInternalConvectionCoef = 5;      % [W/m²/K] : COnvection coefficient between internal air (hive body)
                                  % and surroundings hive walls
lExternalConvectionCoef = 20;     % [W/m²/K] : COnvection coefficient between internal air (hive body)
                                  % and surroundings hive walls
% Vector containing thermal emissivity to sun radiation
lExternalEmissivity = [repmat(lMatwood.eps,1,5) lMatsteal.eps];
lWallLength = 0.02;               % [m] : Hive wall thickness
lNodeCount = 10;                  % Number of wall nodes along the length-axis.

lRoofThickness = 0.0006;          % [m] : Thickness of the roof made of steal
lRoofContactResistance = 0; %0.5E-3/lMatair.lambda;
                                  % [Km²/W] : Contact resistance between steal roof and overframe

lInitTemperature = 25;            % [°C or K] : Initial temperature

tSpec = {0 24*3600/1000 1000};    % [s] : Time vector specification {t0, tf, dt}

lGeoPosition = [43 36 39; ...     % Longitude [deg minute second]
                 3 52 38];        % Latitude  [deg minute second]
lDate = [2023 01 10];             % Date [year, month, day]

% View parameters
lOnlyOneFigure = true;
% Display colors used to plot temperature for each face
% It allows faces to have always the same color
lFaceColorMat = [  0.00000   0.44700   0.74100;... XM
                   0.85000   0.32500   0.09800;... XP
                   0.92900   0.69400   0.12500;... YM
                   0.49400   0.18400   0.55600;... YP
                   0.46600   0.67400   0.18800;... ZM
                   0.30100   0.74500   0.93300;... ZP
                   0.63500   0.07800   0.184000]; % Roof

% Create the figure(1) that will contain the geometry
if lOnlyOneFigure, figure(1); clf; subplot(2,3,1); ...
else figure(1); clf; endif
set(gca, 'FontSize', 14);
xlabel('X', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Y', 'fontsize', 16, 'fontweight', 'bold');
zlabel('Z', 'fontsize', 16, 'fontweight', 'bold');
axis equal;

lDisplayOptions = struct( 'explode', 1.2, ...
                          'explodeCenter', 0.5*lBodySize, ...
                          'displaygrid', true, ...
                          'displayedge', false, ...
                          'facecolortype', 'flat');

% Create model main object with name "main"
lModel = HT_Model_Init('', 'main');

% Create the hive body with dimension 0.5x0.5x0.5m
[lModel0D lModel0D_faces] = HT_Model_Conduction0D("model_hiveBody", ...
                  struct('dim', lBodySize, ...              % Hive body dimension
                         'material', lMatair, ...           % Hive body is filled with air
                         'axis', lGlobalAxis, ...
                         'globalPosition', lGlobalPosition), ...
                  setfield(lOptions, 'nodeName', 'nbody')); % Add the field "nodeName" to the default structure "lOptions"
                                                            % Specify the name of the node (instead of default name)
HT_Plot_Face(lModel0D_faces, lDisplayOptions); % Display the geometry

% Merge the 0D Model to the main model
lModel = HT_Model_Merge(lModel, lModel0D, lOptions);

## Create the 5 walls (except the roof side), XM, XP, YM, YP, ZM
## lSideWallsParams is an array of struct specifying parameter for each walls
## "dim", "axis", "globalposition" can be explicitely specified, but it is easier
## to use field "base". It consists in building a wall by expanding a face.
##                'dim', {  [lBodySize(Y)   lBodySize(Z)], ... % Struct array, XM
##                          [lBodySize(Y)   lBodySize(Z)], ... % XP
##                          [lBodySize(X)   lBodySize(Z)], ... % YM
##                          [lBodySize(X)   lBodySize(Z)], ... % YP
##                          [lBodySize(Y)   lBodySize(X)] }, ... % ZM
##                 'axis', { [-1 0 0; 0 -1 0 ; 0 0 1]', ... % local axis of xm wall
##                          lGlobalAxis, ...
##                          [0 -1 0; 1 0 0; 0 0 1]', ...
##                          [0 +1 0; 1 0 0; 0 0 1]', ...
##                          [0 0 -1; 0 1 0; 1 0 0]' }, ...
##                'globalPosition', {  lGlobalPosition + lBodySize(Y)*lGlobalAxis(:,Y), ...
##                                     lGlobalPosition + lBodySize(X)*lGlobalAxis(:,X), ...
##                                     lGlobalPosition, ...
##                                     lGlobalPosition + lBodySize(Y)*lGlobalAxis(:,Y), ...
##                                     lGlobalPosition }, ...
lSideWallsParams = struct(
                'name', { 'model_hiveWall_xm', ...
                          'model_hiveWall_xp', ...
                          'model_hiveWall_ym', ...
                          'model_hiveWall_yp', ...
                          'model_hiveWall_zm', ...
                          'model_hiveWall_zp'}, ...
                'boundaryname', { {'nwallxm_int', 'nwallxm_out'}, ...
                                  {'nwallxp_int', 'nwallxp_out'}, ...
                                  {'nwallym_int', 'nwallym_out'}, ...
                                  {'nwallyp_int', 'nwallyp_out'}, ...
                                  {'nwallzm_int', 'nwallzm_out'}, ...
                                  {'nwallzp_int', 'nwallzp_out'} });
lSideConnectionParams = struct(
                'name', { 'model_hiveWallConnection_xm', ...
                          'model_hiveWallConnection_xp', ...
                          'model_hiveWallConnection_ym', ...
                          'model_hiveWallConnection_yp', ...
                          'model_hiveWallConnection_zm', ...
                          'model_hiveWallConnection_zp'});

lFaceAirExt = struct(
                'name', { 'face_AirExt_xm', ...
                          'face_AirExt_xp', ...
                          'face_AirExt_ym', ...
                          'face_AirExt_yp', ...
                          'face_AirExt_zm', ...
                          'face_AirExt_zp'}, ...
                'displayOffset', {[-0.1 0 0], ...
                                  [+0.1 0 0],...
                                  [0   -0.1  0],...
                                  [0   +0.1  0],...
                                  [0   0    -0.1], ...
                                  [0   0    +0.1]} );

lAirExtConnectionModel = struct(
                'name', { 'model_airExtConnection_xm', ...
                          'model_airExtConnection_xp', ...
                          'model_airExtConnection_ym', ...
                          'model_airExtConnection_yp', ...
                          'model_airExtConnection_zm', ...
                          'model_airExtConnection_zp'});

lWallFaces = cell(6,1); % Store the faces of each 5 walls around the hive
lWallNodes = cell(6,1); % Each cell will contain the list of nodes of the corresponding wall
                        % xm, xp, ym, yp zm

for k=1:6
  p = lSideWallsParams(k);
  % Create a 1-dimensionnal single layer wall model in the XP (X+) direction
  [lModelWall1D lModelWall1D_faces lModelWall1D_nodes] = HT_Model_Conduction1D(...
              p.name, ...
              'length',    lWallLength,             ... % Thickness of the wall
              'material',  lMatwood,                ... % Material
              'base',      lModel0D_faces(k),       ... % Base face
              'n',         lNodeCount,              ... % Spatial discretization. Wall is divided into "n" nodes
              'gridType',   'hh', ...
              'boundaryNames', { p.boundaryname },  ...
              lOptions);
  lWallNodes{k} = lModelWall1D_nodes;
  lWallFaces{k} = lModelWall1D_faces;
  % Merge the 1D Model to the main model
  lModel = HT_Model_Merge(lModel, lModelWall1D, lOptions);

  % Display the 6 faces of the 1D Model
  HT_Plot_Face(lModelWall1D_faces, lDisplayOptions);

  % Create a thermal connection between hive body and surrounding wall
  % Connection with wall XP
  lModel_BodyWall = HT_Model_ConnectFaces(...
                          lSideConnectionParams(k).name,        ... % Model name
                          lModel0D_faces(k),                    ... % First face is the hive body air cube...
                          lModelWall1D_faces(FaceXM),           ... % Second face is the internal face of the wall (here, XM is always the internal wall)
                          struct('g', lInternalConvectionCoef), ... % Internal convection coefficient
                          setfield(lOptions, "checkGeometry", lCheckGeometry)); % Add a field to the structure lOptions
  % Merge the model of the thermal connection to the main model
  lModel = HT_Model_Merge(lModel, lModel_BodyWall, lOptions);
endfor

% Add the roof over the wood face ZP
[lModelRoof1D lModelRoof1D_faces lModelRoof1D_nodes] = HT_Model_Conduction1D(...
                  'model_hiveRoof',                       ...
                  'length',       lRoofThickness,                  ... % Thickness of the wall
                  'material',     lMatsteal,                       ... % Material
                  'base',         lWallFaces{FaceZP}(FaceXP),      ... %
                  'n',            1,                               ... % Spatial discretization. Wall is divided into "n" nodes
                  'boundaryNames', {{ 'nroof_out', 'nroof_in'}},      ...
                  lOptions);
% Merge the roof Model to the main model
lModel = HT_Model_Merge(lModel, lModelRoof1D, lOptions);
clear lModelRoof1D;

% Create a thermal connection between steal roof and overframe
lModel_HiveRoof = HT_Model_ConnectFaces(...
                        'model_hiveWallRoofConnection', ...
                        lModelRoof1D_faces(FaceXM),             ... % First face is the roof internal face
                        lWallFaces{FaceZP}(FaceXP),             ... % Second face is the internal face of the wall (here, XM is always the internal wall)
                        struct('g', 1/lRoofContactResistance),  ... % Internal convection coefficient
                        setfield(lOptions, "checkGeometry", lCheckGeometry)); % Add a field to the structure lOptions
% Merge the model of the thermal connection to the main model
lModel = HT_Model_Merge(lModel, lModel_HiveRoof, lOptions);


% Display the 6 faces of the 1D Model
HT_Plot_Face(lModelRoof1D_faces, setfield(lDisplayOptions, 'explodeOffset', [0.0 0.0 0.05]));

% Creation of a node for external air
lAirExtNode = HT_Node_Init('AirExt');
% Connect XM,XP,YM,YP,ZM (not ZP) to outside air and Roof.ZP to outside air
lExternalFaces = cellfun(@(v) v(FaceXP), lWallFaces(1:5), 'UniformOutput', false); % Retrieve the external faces from cell <lWallFaces>
% Add the external face of the roof to the list of external faces
lExternalFaces = [lExternalFaces; lModelRoof1D_faces(FaceXP)];
% Connect external faces to the outside air
lNodeModel = HT_Model_Connect('Model_AirExtConnection',       ...
                              lAirExtNode,                    ...
                              lExternalFaces,                 ... Cell array which contains the external faces of the hive
                              'g', lExternalConvectionCoef,   ... Uniform conductance between external faces and outside air
                              lOptions);
lModel = HT_Model_Merge(lModel, lNodeModel, lOptions);
clear lExternalFaces lNodeModel lAirExtNode ;



% Creation of a node for sky temperature
##lSkyNode = HT_Node_Init('SkyExt');
##lNodeModel = HT_Model_Connect('Model_SkyConnection', ...
##                              lSkyNode, ...
##                              cellfun(@(v) v(FaceXP), lWallFaces, 'UniformOutput', false), ... % Retrieve the external faces from cell <lWallFaces>
##                              'g', {2 2 2 2 0 4}, ... Uniform conductance between external faces and outside air
##                              lOptions);
##lModel = HT_Model_Merge(lModel, lNodeModel);


clear lModelWall1D lModelWall1D_nodes lModel_BodyWall p lModel_WallAirExt;
clear lSideWallsParams lSideConnectionParams lFaceAirExt lAirExtConnectionModel;


% ======================= Command definition ==========================
% Retrieve the sun position at current date time
[zenithAngle azimuthAngle] = HT_GetSolarAngle(lDate(1),                           ... Year
                                              lDate(2:3),                         ... [Month Day]
                                              tSpec{1} + (0:(tSpec{3}-1))*tSpec{2}, ... Day time
                                              lGeoPosition(1,:)',                 ... Longitude
                                              lGeoPosition(2,:)');                %  Latitude

% Retrieve the cosine of angle between each hive faces and the sun
[cosBeta] = HT_GetHiveSolarAngles(...
                          zenithAngle,  ... % Sun zenith angle (90°=sun at horizon)
                          azimuthAngle, ... % Sun azimuth angle (radian)
                          [lHiveAzimuth,lHiveTilt], ...
                          struct('angle', true)); % angle=true means previous
                               %  argument is expressed as [azimuth tilt]

% Retrieve the exposition coefficient (negative values are clamped to 0)
lExpositionMat = HT_GetSolarFluxFactor(cosBeta, lExternalEmissivity);

% Nullify exposition coefficients when the sun is below the horizon
% i.e. when zenith angle is higher than 90°=pi/2
lExpositionMat(zenithAngle >= pi/2, :) = 0;
lZenithAttenuation = 0.8; % Voir Solar Irradiance on Wikipedia (radiation at ground level over/extra atmosphere) with at zenith
lSolarIrradiance = 1361 * exp(log(lZenithAttenuation) ./ max(0.0, cos(zenithAngle)));
% Multiply each line by the solar irradiance
lExpositionMat = lSolarIrradiance .* lExpositionMat;
% Retrieve the area of each external faces
lAreaFaces = HT_Face_GetAbsoluteData(cellfun(@(v) v(FaceXP), lWallFaces, 'UniformOutput', false), 'nodesArea');
% Multiply each columns by the face area
lSolarRadiationMat = lExpositionMat .* cell2mat(lAreaFaces)';

lCmd = HT_Cmd_Init(
        'name', 'cmd_airext', ...
          'type', 'temperature', ...
          'nodes', 'AirExt', ...
          'data', 25, ... External temperature of 10°C
##        'name', 'cmd_skyext', ...
##          'type', 'temperature', ...
##          'nodes', 'SkyExt', ...
##          'data', -50, ... External temperature of 10°C
        'name', 'cmd_sun_xm', ...
          'type', 'flux', ...
          'nodes', 'nwallxm_out', ...
          'data', lSolarRadiationMat(:,FaceXM)', ...
        'name', 'cmd_sun_xp', ...
          'type', 'flux', ...
          'nodes', 'nwallxp_out', ...
          'data', lSolarRadiationMat(:,FaceXP)', ...
        'name', 'cmd_sun_ym', ...
          'type', 'flux', ...
          'nodes', 'nwallym_out', ...
          'data', lSolarRadiationMat(:,FaceYM)', ...
        'name', 'cmd_sun_yp', ...
          'type', 'flux', ...
          'nodes', 'nwallyp_out', ...
          'data', lSolarRadiationMat(:,FaceYP)', ...
        'name', 'cmd_sun_zp',     ...
          'type', 'flux',         ...
          'nodes', 'nroof_out', ... Roof face
          'data', lSolarRadiationMat(:,FaceZP)'...
          ); ... Solar irradiance of 1000W/m²

% ============================== Model resolution ==============================
[T, U, Tnodes, Unodes, t] = HT_SolveModel(  ...
         lModel,                            ...  Model to be solved
         lCmd,                              ...    List of commands
         lInitTemperature,                  ... Initial temperature
         tSpec,                             ... Time vector
         struct('all', true,                ...
                'algorithm', 'chol',        ...
                'replaceT0', true,          ...
                'verbose', lOptions.verbose,...
                'info', 'default'));

% ======================= Display temperature evolution ========================
if lOnlyOneFigure, subplot(2,3,[2 3]); else figure(2); clf; endif
% Extract the temperature evolution of node called "nbody"
tplot = t(1) + (0:0.01:1)*(t(end)-t(1));
lPlotInfos = {'AirExt'      , 'Air ext.',       '-',  'black'                 , 4; ...
              'nbody'       , 'Hive body',      '-' , 'red'                   , 4; ... List of nodes whose temperature
              'nwallxm_int' , 'XM face (int.)', '--', lFaceColorMat(FaceXM,:) , 2; ... has to be retrieved and
              'nwallxm_out' , 'XM face (ext.)', '-' , lFaceColorMat(FaceXM,:) , 2; ... corresponding plot name
              'nwallxp_int' , 'XP face (int.)', '--', lFaceColorMat(FaceXP,:) , 2; ...
              'nwallxp_out' , 'XM face (ext.)', '-' , lFaceColorMat(FaceXP,:) , 2; ...
              'nwallym_int' , 'YM face (int.)', '--', lFaceColorMat(FaceYM,:) , 2; ...
              'nwallym_out' , 'YM face (ext.)', '-' , lFaceColorMat(FaceYM,:) , 2; ...
              'nwallyp_int' , 'YP face (int.)', '--', lFaceColorMat(FaceYP,:) , 2; ...
              'nwallyp_out' , 'YP face (ext.)', '-' , lFaceColorMat(FaceYP,:) , 2; ...
              'nwallzm_int' , 'ZM face (int.)', '--', lFaceColorMat(FaceZM,:) , 2; ...
              'nwallzm_out' , 'ZM face (ext.)', '-' , lFaceColorMat(FaceZM,:) , 2; ...
              'nwallzp_int' , 'ZP face (int.)', '--', lFaceColorMat(FaceZP,:) , 2; ...
              'nwallzp_out' , 'ZP face (ext.)', '-' , lFaceColorMat(FaceZP,:) , 2; ...
##              'nroof_int'   , 'Roof (int.)'   , '--', lFaceColorMat(7,:)      , 2; ... Only used if more than 1 node is defined in the roof
              'nroof_out'   , 'Roof (ext.)'   , '-' , lFaceColorMat(7,:)      , 2};
lPlotTemp = HT_Result_GetNodeTemperature(...
                    lPlotInfos(:,1),  ... List of nodes
                    T,                ... Temperature vector
                    Tnodes,           ... Model modes list
                    'time', tplot,    ... Display time vector
                    'timevector', t); ... Time Vector
for i=1:(rows(lPlotInfos))
  plot(tplot/3600,  lPlotTemp(i,:),               ... Temperature data
                    lPlotInfos{i,3},              ... Line style
                    'linewidth', lPlotInfos{i,5}, ... Line width
                    'color', lPlotInfos{i,4},     ... Line color
                    'displayname', lPlotInfos{i,2}); hold on;
endfor
set(gca, 'fontsize', 14);
xlabel('Temps [hour]');
ylabel('Temperature [degC]');
title('Temperature of hive body and of int./ext. face of each side');
hl = legend();
set(hl, 'fontsize', 10);
grid on;

% ================ Display temperature profils inside Xp wall ==================
if lOnlyOneFigure, subplot(2,3,4); else figure(3); clf; endif
% Extract temperatures for each nodes of XP wall. The list of nodes was saved
% in <lWallNodes{FaceXP}>. The time is interpolated and is regularly spaced.
##tprofil = logspace(0, log10(t(end)), 10); %0:600:t(end); % [s]
tprofil = 0:7200:t(end); % [s]
nNodes = numel(lWallNodes{FaceXP});
Twallprofil = HT_Result_GetNodeTemperature(lWallNodes{FaceXP}, T, Tnodes, ...
                                          'time', tprofil, ... Time of different profils
                                          'timevector', t);
for k=1:size(Twallprofil,2)
  h = plot(Twallprofil(:,k), 'linewidth', 2.0); hold on;
  text(nNodes/2, Twallprofil(idivide(nNodes,int32(2)),k), ...
        sprintf('t=%.1fh', tprofil(k)/3600),...
        'fontsize', 14, ...
        'fontweight', 'bold');
endfor

set(gca, 'fontsize', 14);
xlabel('Position [node]');
ylabel('Temperature [degC]');
title('Temperature profils inside hive wall XP (north)');
grid on;

% ================ Display temperature profils inside Xm wall ==================
if lOnlyOneFigure, subplot(2,3,5); else figure(4); clf; endif
% Extract temperatures for each nodes of XP wall. The list of nodes was saved
% in <lWallNodes{FaceXP}>. The time is interpolated and is regularly spaced.
##tprofil = logspace(0, log10(t(end)), 10); %0:600:t(end); % [s]
tprofil = 0:7200:t(end); % [s]
nNodes = numel(lWallNodes{FaceXM});
Twallprofil = HT_Result_GetNodeTemperature(lWallNodes{FaceXM}, T, Tnodes, ...
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
title('Temperature profils inside hive wall XM (south)');
grid on;

% =================== Display solar irradiance for each face ===================
if lOnlyOneFigure, subplot(2,3,6); else figure(5); clf; endif
hold on;
lFaceName = { 'XM', 'XP', 'YM', 'YP', 'ZM', 'ZP' };
for i=1:6
  h = plot(t/3600,  lExpositionMat(:,i),              ...
                    'linewidth', 3.0,                 ...
                    'color', lFaceColorMat(i,:),  ...
                    'displayname', lFaceName{i});
  endfor
plot(t/3600, lSolarIrradiance, '--k', 'linewidth', 3.0, 'displayname', 'Solar irradiance');
ylabel('Solar irradiance [W/m²]', 'Fontsize', 14);
xlabel('Time [hour]', 'Fontsize', 14);
legend();
title('Direct solar irradiance incoming on each hive faces', 'Fontsize', 14);
grid on;
