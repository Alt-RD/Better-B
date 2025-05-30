%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University, AltRD-Emmanuel Ruffio
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
% Script d'initialisation des variables globales utilisées par le logiciel HiveTemp

X=1; Y=2; Z=3;
DIRstr = {'X', 'Y', 'Z'};

U=1; V=2;

FaceXM=1; FaceXP=2; FaceYM=3; FaceYP=4; FaceZM=5; FaceZP=6;
FaceDIRstr = {'XM', 'XP', 'YM', 'YP', 'ZM', 'ZP'};
FaceDIRstrL = {'xm', 'xp', 'ym', 'yp', 'zm', 'zp'};
FaceDIRcoef = {'1', '-1', '1', '-1', '1', '-1'};
FaceDIRmat = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1]';

HT_WAVELENGTH_VISIBLE = 1;    % Around 1um
HT_WAVELENGTH_IR_LOW = 2;
HT_WAVELENGTH_IR_MED = 10;    % Around 10um

% Boundary condition type
HT_BTYPE_T_CONSTANT_UNIFORM = 11;
HT_BTYPE_T_CONSTANT_NON_UNIFORM = 12;
HT_BTYPE_T_VARIABLE_UNIFORM = 13;
HT_BTYPE_T_VARIABLE_NON_UNIFORM = 14;
HT_BTYPE_T_INTERPOLATION_ARRAY = 15;

HT_BTYPE_FLUX_CONSTANT_UNIFORM = 21;
HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM = 22;
HT_BTYPE_FLUX_VARIABLE_UNIFORM = 23;
HT_BTYPE_FLUX_VARIABLE_NON_UNIFORM = 24;
HT_BTYPE_FLUX_INTERPOLATION_ARRAY = 25;

global HT_VAR_EPSILON_POS;
global HT_VAR_EPSILON_U;
global HT_VAR_CHECK_INPUT;
global HT_VAR_STRICT_INPUT_PARAMETERS;



