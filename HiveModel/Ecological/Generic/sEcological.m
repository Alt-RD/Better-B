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
% See Model_Ecological_Schemas.pdf
%
% Todo and improvements:
% 1) Improve geometry between each sides
% 2) Add thermal connexion between each sides
% 3) Add radiation heat transfer inside the hive and in the overframe
% 4) Improve geometry above the overframe (the wall thickness changes)
% ========================================================================

clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../../../Better-B/HiveTemp/source/';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

% ========================================================================
%                              PARAMETERS
% ========================================================================
sEcological_setupParams();

% ========================================================================
%                            HIVE THERMAL MODEL
% ========================================================================
sEcological_setupModel();

% ========================================================================
%                               COMMANDS
% ========================================================================
sEcological_setupCommands();

% ========================================================================
%                               RESOLUTION
% ========================================================================
sEcological_solve();

##sTest

