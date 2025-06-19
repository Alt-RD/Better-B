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
% This is a program to test the function "Standalone_Solve_Hygro1D".
% The water diffusion equation is solved to get the water content of a wood
% sample during a drying process.

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

lExtConditions = struct(      'name', {'left', 'right'}, ...
                              'T', {20, 30}, ...
                              'RH', {0.6, 0.6});

lSorptionf = @(RH, T) 0.244*exp(0.69*log(RH).*exp(1.69*RH));
lWeqVector = lSorptionf([lExtConditions.RH], [lExtConditions.T]);

lParamsNLCase = struct(...
                 'length', 0.05,                      ... [cm] length of the wall
                 'diffusion', @(t,i,W) 2.1E-11*exp(5.77*W),   ... [function_handle D = @(W)] diffusion coefficient
                 'convection', [1E-5, 1E-5],          ... [1x2] or [nt x 2] convection coefficient for each sides of the wall
                 'weq', lWeqVector,                   ... [1x2] or [nt x 2] equilibrium water content on each sides
                 'winit', 0.1,                        ... [1x1] or [nx1] initial water content
                 'n', 21, ...                         ... number of nodes
                 'gridType', 'ff', ...
                 't', (0:0.01:1)'.^2*200*24*3600);

% For the linear case, the diffusion coefficient is constant
lParamsLinearCase = setfield(lParamsNLCase, 'diffusion', @(t,i,W) repmat(2.1E-11*exp(5.77*0.1), size(W)));

lOptions = struct('verbose', true);
[WMatNL lInfosNL] = HT_Standalone_Solve_Hygro1D(lParamsNLCase, lOptions);

[WMatL lInfosL] = HT_Standalone_Solve_Hygro1D(lParamsLinearCase, lOptions);

sTest_StandaloneSolveHygro1D_Plot;

