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



% ========================================================================
%                          SOLVE THE MODEL
% ========================================================================

% Time vector

lStartTimeSec = sum(lHiveParams.time .* [3600 60 1]);

[T, U, Tnodes] = HT_SolveModel(...
         lModel,                                            ... Model to be solved
         lCmd,                                              ... List of commands
         lComputation.initTemperature,                      ... Initial temperature
         { lStartTimeSec, lComputation.timeStep, lComputation.nt},  ... Time vector
         struct('all', true,                                ... All temperature must be returned. T will be (dim Nnodes x Ntimes)
                'algorithm', 'linear',                      ... Algorithm used
                'replaceT0', true,                          ... Overwrite all initial temperature with specified temperature <lInitTemperature>
                'verbose', lOptions.verbose,                ...
                'unit', 'degres',                           ... % Warns the function that all temperature are expressed in degC.
                'info', 'default'));

##sTest;

##figure(2);
##clf;
##
##HT_Plot_Face(lMod_CylinderFaces.outside, 'temperature', T(:,end), 'nodes', Tnodes, l3DPlotOptions);
##HT_Plot_Face(lMod_CylinderFaces.inside, 'temperature', T(:,end), 'nodes', Tnodes, l3DPlotOptions);




