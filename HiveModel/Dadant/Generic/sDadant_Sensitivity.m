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
% Load hive parameters
sDadant_setupParams;

% Load parameters related to sensitivity study
% The file containing the parameters is specified by "lSensitivityParams.parameterFile"
lSetupFile = lSensitivityParams.parameterFile;
if endsWith(lSetupFile, '.m'), lSetupFile((end-1):end) = []; endif;
eval(strcat(lSetupFile, ';'));

lCacheReset = true; % If true, the cache file is erased.

% ========================================================================
%                          SENSITIVITY
% ========================================================================
nParam = numel(lSensitivityAnalysis);

lParameterData = struct(...
  'name', {lSensitivityAnalysis.name});

% Sensitivity matrix
XInfos = struct('name', '', ...
                'X', repmat({{}}, nParam, 1),
                'T', [], ...
                't', [], ...
                'faces', [], ...
                'volumes', [], ...
                'nodes', []);

lAbsoluteDependencyFile = cellfun(@(v) make_absolute_filename(v), lSensitivityParams.dependencies, 'UniformOutput', false);
lFileStats = struct('file', lSensitivityParams.dependencies, ...
                    'absoluteFile', lAbsoluteDependencyFile, ...
                    'stat', cellfun(@(v) stat(v), lAbsoluteDependencyFile, 'UniformOutput', false));

do
  % Attempt to reload data from cache file
  lCacheAbsoluteFilePath = lSensitivityParams.cacheFile;

  if ~isempty(lCacheAbsoluteFilePath)
    lCacheAbsoluteFilePath = make_absolute_filename(lCacheAbsoluteFilePath);
  endif
  if ~lCacheReset && ~isempty(lCacheAbsoluteFilePath) && isfile(lCacheAbsoluteFilePath)
    lCacheData = load(lCacheAbsoluteFilePath);
    if ~isfield(lCacheData, 'fileStats'), break; endif
    if ~isequal(lCacheData.fileStats, lFileStats), break; endif
    if ~isfield(lCacheData, 'sensitivityParams') || ~isequal(lCacheData.sensitivityParams, lSensitivityParams), break; endif

    for i=1:nParam
      % Try to find the parameter name in the cache data
      ind = find(strcmpi({lCacheData.XInfos.name}, lSensitivityAnalysis(i).name));
      if ~isempty(ind)
        XInfos(i) = lCacheData.XInfos(ind);
        disp(sprintf('=> Reload data from backup file: Parameter <%s>', lParameterData(i).name));
      endif
    endfor
  endif
until true;

for i=1:nParam
  if ~isempty(XInfos(i).name), continue; endif; % Loaded from cache file ?

  disp(sprintf('=> Sensitivity analysis of parameter <%s>', lParameterData(i).name));

  % Compute sensitivity using finite difference
  % Parameters of the "prev"/"next" model with slight change of one parameter
  lParamsPrev = lSensitivityAnalysis(i).setterPrev(lSensitivityAnalysis(i), lParams);
  lParamsNext = lSensitivityAnalysis(i).setterNext(lSensitivityAnalysis(i), lParams);

  lOptions.verbose = false;

  % Build the thermal model
  [lModelP lFaceStructP lVolumeStructP] = Dadant_setupModel(lParamsPrev);
  [lModelN lFaceStructN lVolumeStructN] = Dadant_setupModel(lParamsNext);

  % Solve the thermal model
  t = lComputation.startTime + (0:(lComputation.nt-1)) * lComputation.timeStep;

  % Setup commands
  assert(isfile(strcat(lParamsPrev.command.scriptFile, '.m')), 'The specified file for command is not valid');
  eval(strcat('lCmd = ', lCommand.scriptFile, '(lParamsPrev, lFaceStructP, lOptions);'));

  [TP, ~, TnodesP] = HT_SolveModel(...
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

  [TN, ~, TnodesN] = HT_SolveModel(...
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

  assert(all(strcmp(TnodesN, TnodesP)));

  % Save results in the XInfos structure
  XInfos(i).X = lSensitivityAnalysis(i).compute(lSensitivityAnalysis(i), lParams, TP, TN);
  XInfos(i).name = lParameterData(i).name;
  XInfos(i).t = t;
  XInfos(i).T = (TP+TN)/2;
  XInfos(i).nodes = TnodesN;
  XInfos(i).faces = lFaceStructN;
  XInfos(i).volumes = lVolumeStructN;
endfor

if ~isempty(lCacheAbsoluteFilePath)
  sensitivityParams = lSensitivityParams;
  fileStats = lFileStats;

  save('-binary', lCacheAbsoluteFilePath, 'XInfos','sensitivityParams', 'fileStats');
  clear sensitivityParams fileStats;
endif

sDadant_SensitivityPlot;

