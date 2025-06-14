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
%
% See Model_Dadant_Schemas.pdf
% This script computes the sensitivity of parameters specified below
% ========================================================================

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../../HiveTemp/source/';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

% ========================================================================
%                          USER PARAMETERS
% ========================================================================
sDadant_setupParams;

% ========================================================================
%                          SENSITIVITY
% ========================================================================
##eval(lHiveParamSetupScript);      % Setup Parameters

nParam = numel(lSensitivityAnalysis);

lParameterData = struct(...
  'name', {lSensitivityAnalysis.name});

% Sensitivity matrix
XInfos = struct('name', '', ...
                'X', repmat({{}}, nParam, 1),
                'T', [], ...
                'faces', [], ...
                'volumes', []);

for i=1:nParam
  disp(sprintf('=> Sensitivity analysis of parameter <%s>', lParameterData(i).name));

  % Compute sensitivity using finite difference
  % Parameters of the "prev"/"next" model with slight change of one parameter
  lParamsPrev = lSensitivityAnalysis(i).setterPrev(lParams);
  lParamsNext = lSensitivityAnalysis(i).setterNext(lParams);

  lOptions.verbose = false;

  % Build the thermal model
  [lModelP lFaceStructP lVolumeStructP] = Dadant_setupModel(lParamsPrev);
  [lModelN lFaceStructN lVolumeStructN] = Dadant_setupModel(lParamsNext);

  % Solve the thermal model
  t = lComputation.startTime + (0:(lComputation.nt-1)) * lComputation.timeStep;

  % Setup commands
  assert(isfile(strcat(lParamsPrev.command.scriptFile, '.m')), 'The specified file for command is not valid');
  eval(strcat('lCmd = ', lCommand.scriptFile, '(lParamsPrev, lFaceStructP, lOptions);'));

  [TP] = HT_SolveModel(...
           lModelP,                                            ... Model to be solved
           lCmd,                                              ... List of commands
           lComputation.initTemperature,                      ... Initial temperature
           { lComputation.startTime lComputation.timeStep, lComputation.nt},  ... Time vector
           struct('all', true,                                ... All temperature must be returned. T will be (dim Nnodes x Ntimes)
                  'algorithm', 'linear',                      ... Algorithm used
                  'replaceT0', true,                          ... Overwrite all initial temperature with specified temperature <lInitTemperature>
                  'verbose', lOptions.verbose,                ...
                  'unit', 'degres',                           ... % Warns the function that all temperature are expressed in degC.
                  'info', 'default'));

  assert(isfile(strcat(lParamsNext.command.scriptFile, '.m')), 'The specified file for command is not valid');
  eval(strcat('lCmd = ', lCommand.scriptFile, '(lParamsNext, lFaceStructN, lOptions);'));

  [TN] = HT_SolveModel(...
           lModelN,                                            ... Model to be solved
           lCmd,                                              ... List of commands
           lComputation.initTemperature,                      ... Initial temperature
           { lComputation.startTime lComputation.timeStep, lComputation.nt},  ... Time vector
           struct('all', true,                                ... All temperature must be returned. T will be (dim Nnodes x Ntimes)
                  'algorithm', 'linear',                      ... Algorithm used
                  'replaceT0', true,                          ... Overwrite all initial temperature with specified temperature <lInitTemperature>
                  'verbose', lOptions.verbose,                ...
                  'unit', 'degres',                           ... % Warns the function that all temperature are expressed in degC.
                  'info', 'default'));

  lUserData = lSensitivityAnalysis(i).getUserData(lParams);

  % Save results in the XInfos structure
  XInfos(i).X = lSensitivityAnalysis(i).compute(lUserData, TP, TN);
  XInfos(i).name = lParameterData(i).name;
  XInfos(i).T = (TP+TN)/2;
  XInfos(i).nodes = lModelN.nodes;
  XInfos(i).faces = lFaceStructN;
  XInfos(i).volumes = lVolumeStructN;
endfor

sDadant_SensitivityPlot;

