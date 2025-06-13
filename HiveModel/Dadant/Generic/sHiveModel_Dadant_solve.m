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
%  This work was supported by the Better-B project, which has received
%  funding from the European Union, the Swiss State Secretariat for
%  Education, Research and Innovation (SERI) and UK Research and
%  Innovation (UKRI) under the UK government's Horizon Europe funding
%  guarantee (grant number 10068544).'
%
% See Model_Dadant_Schemas.pdf
% ========================================================================

assert(exist('lCmd', 'var') == 1, 'No command structure found. Command script must be executed first');

% ========================================================================
%                               RESOLUTION
% ========================================================================

% Time vector
t = lComputation.startTime + (0:(lComputation.nt-1)) * lComputation.timeStep;

[T, U, Tnodes] = HT_SolveModel(...
         lModel,                                            ... Model to be solved
         lCmd,                                              ... List of commands
         lComputation.initTemperature,                      ... Initial temperature
         { lComputation.startTime lComputation.timeStep, lComputation.nt},  ... Time vector
         struct('all', true,                                ... All temperature must be returned. T will be (dim Nnodes x Ntimes)
                'algorithm', 'linear',                      ... Algorithm used
                'replaceT0', true,                          ... Overwrite all initial temperature with specified temperature <lInitTemperature>
                'verbose', lOptions.verbose,                ...
                'unit', 'degres',                           ... % Warns the function that all temperature are expressed in degC.
                'info', 'default'));


% ============================== Model resolution ==============================




