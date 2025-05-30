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
% This script shows the heat radiation inside a empty box. The radiation
% luminance matrix is retrieved from function HT_Model_RadiationInBox.
% No model is solved. The radiation matrix is just displayed to shows how
% a node radiates in all direction.
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
lMatwood = HT_Material_New('wood', 'eps',    [HT_WAVELENGTH_VISIBLE  0.5; ...
                                              HT_WAVELENGTH_IR_MED   0.5]);

% Defines the global axis (orientation of the hive) and position
lGlobalAxis = eye(3);              % Identity matrix denotes the 3 main direction of X,Y,Z in columns
lGlobalPosition = [0.0; 0.0; 0.0]; % [m x m x m] Position of the hive

% Global parameters
lOptions = struct('verbose', true);  % Display additional information (level 4)
lCheckGeometry = true;            % Perform additional tests to check geometry coherency
                                  % (in HT_Model_ConnectFaces function)

lBodySize = [0.5 0.5 0.5]';       % [m x m x m]: Hive body size along X,Y,Z axis direction

lInternalConvectionCoef = 5;      % [W/m²/K] : COnvection coefficient between internal air (hive body)
                                  % and surroundings hive walls

lNodeCount = int32(7);                  % Number of wall nodes along the length-axis.

lInitTemperature = 25;            % [°C or K] : Initial temperature

tSpec = {0 24*3600/1000 1000};    % [s] : Time vector specification {t0, tf, dt}

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

% Create model main object with name "main"
lModel = HT_Model_Init('', 'main');

% Create the six faces of the box, XM,XP,YM,YP,ZM,ZP
lBoxParams = struct('name', { 'fxm', 'fxp', 'fym', 'fyp', 'fzm', 'fzp' }, ...
                    'size', {   lBodySize([Y Z]), ...
                                lBodySize([Y Z]), ...
                                lBodySize([X Z]), ...
                                lBodySize([X Z]), ...
                                lBodySize([X Y]), ...
                                lBodySize([X Y]) }, ...
                    'axis', {   [-1 1] .* lGlobalAxis(:, [Y Z]), ...
                                [+1 1] .* lGlobalAxis(:, [Y Z]), ...
                                [+1 1] .* lGlobalAxis(:, [X Z]), ...
                                [-1 1] .* lGlobalAxis(:, [X Z]), ...
                                [-1 1] .* lGlobalAxis(:, [X Y]), ...
                                [+1 1] .* lGlobalAxis(:, [X Y]) }, ...
                    'norm', {   -lGlobalAxis(:, X), ...
                                +lGlobalAxis(:, X), ...
                                -lGlobalAxis(:, Y), ...
                                +lGlobalAxis(:, Y), ...
                                -lGlobalAxis(:, Z), ...
                                +lGlobalAxis(:, Z) }, ...
                    'globalPosition', { [ 0 1 0 ]' .* lBodySize, ...
                                        [ 1 0 0 ]' .* lBodySize, ...
                                        [ 0 0 0 ]' .* lBodySize, ...
                                        [ 1 1 0 ]' .* lBodySize, ...
                                        [ 1 0 0 ]' .* lBodySize, ...
                                        [ 0 0 1 ]' .* lBodySize }, ...
                    'material', { lMatwood, ...
                                  lMatwood, ...
                                  lMatwood, ...
                                  lMatwood, ...
                                  lMatwood, ...
                                  lMatwood } );

% Create the six faces
lBoxFaces = HT_Face_Init(lBoxParams);

% Create a regular mesh (lNodeCount x lNodeCount) for each faces
lBoxFaces = HT_Face_CreateMesh(lBoxFaces,   'n', lNodeCount, ...
                                            'type', 'uniform', ...
                                            'uEdgeType', 'ff', ...
                                            'vEdgeType', 'ff');

% Create the radiation model and retrieve the fluxMatrix [-]
% For a given vector of temperature T, the amount of incoming radiation
% for each nodes is given by Phi[W/m²] = lFluxMatrix * sigma * T^4
% HT_SolveModel does the same but multiply by node areas to get the absolute
% heat flux [W]
[lModel_Radiation lFluxMatrix] = HT_Model_RadiationInBox('mRadBox', ...
                                            struct( 'face', lBoxFaces, ...
                                                    'wavelength', HT_WAVELENGTH_IR_MED));

figure(1);
clf;

lNodeIndex = { 1 2, ...
               lNodeCount^2*5+idivide(lNodeCount, 2)*lNodeCount+idivide(lNodeCount, 2)+1 , ...
               lNodeCount^2*5+(idivide(lNodeCount, 2))*lNodeCount+1 };

for i=1:4
  subplot(2,4,i);

  T = repmat(273.15, rows(lFluxMatrix), 1);
  T(lNodeIndex{i}) += 1;

  lFluxVector = lFluxMatrix * (5.67E-8 * T.^4);
  lFluxVector /= (5.67E-8 * (274.15^4-273.15^4));
  lFluxVector(lNodeIndex{i})= NaN;
##  lFluxVector = lFluxMatrix(:,lNodeIndex{i});
##  lFluxVector(lNodeIndex{i})= NaN;

  campos([-3.2809   2.2290   1.9908]);
  camup([0 0 1]);
  HT_Plot_Face(lBoxFaces([1, 4, 6]), 'temperature', log10(lFluxVector), 'nodes', lModel_Radiation.rad{1}.nodes, l3DPlotOptions);
  set(gca, 'FontSize', 14);
  xlabel('x');
  ylabel('y');
  zlabel('z');
  axis('equal');
  ##colorbar();
  title(sprintf('Intensity [dB] (top) (\\epsilon=%.1f)', HT_Material_GetEmissivity(lMatwood, HT_WAVELENGTH_IR_MED)), 'fontsize', 16, 'fontweight', 'bold');
  grid on;
  caxis([-4 -0.5]);
  colorbar();

  % =======================================
  subplot(2,4,4+i);

  campos([-3.3311   2.4228  -1.7179]);
  camup([0 0 1]);
  HT_Plot_Face(lBoxFaces([1, 4, 5]), 'temperature', log10(lFluxVector), 'nodes', lModel_Radiation.rad{1}.nodes, l3DPlotOptions);
  set(gca, 'FontSize', 14);
  xlabel('x');
  ylabel('y');
  zlabel('z');
  axis('equal');
  ##colorbar();
  title('(bottom)', 'fontsize', 16, 'fontweight', 'bold');
  grid on;
  caxis([-4 -0.5]);
  colorbar();
endfor

