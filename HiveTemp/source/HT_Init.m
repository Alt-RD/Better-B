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
%> @file  HT_Init.m
%> @brief Initialize the HiveTemp library and import constant symbols
%>
%> This function performs severals initialization necessary for the
%> library to work. The path variable of Octave is especially updated
%> to include all directories containing HiveTemp functions. The octave
%> variable <HT_VAR_LIB_PATH> must be set before calling this function.
%> However, since HT_Init() function must be reachable by Octave, it is
%> necessary to add manually the main directory of HiveTemp to the path
%> variable of Octave. As a result, the standard way to create a new model
%> with HiveTemp is as follow:
%>
%>@code
%>clear variables;
%>% Specify the main directory of HiveTemp library (relative or absolute)
%>HT_VAR_LIB_PATH = './../HiveTemp/';
%>addpath(make_absolute_filename(HT_VAR_LIB_PATH));
%>
%>% Init the library
%>HT_Init();
%>@endcode
%>
%> This file contains global variables used by the library. This variable
%> can be further modified in the user script. See @ref HT_VAR_EPSILON_POS
% ======================================================================

%======================================================================
%> @mainpage Documentation for HiveTemp library
%>
%> @section intro Introduction
%>
%> The @b HiveTemp library (https://www.alt-rd.com) allows you to create thermal
%> modeling of simple objects.
%>
%> See @ref HT_Init.m file to create a new model.
%>
%>
%>@verbatim
%>% ======================================================================
%>%> @brief Brief description of the function
%>%>
%>%> More detailed description.
%>%>
%>%> @param arg1 First argument
%>%> @param arg2 Second argument
%>%>
%>%> @retval out1 return value for the first output variable
%>%> @retval out2 return value for the second output variable
%>% ======================================================================
%>[out1, out2] = function( arg1, arg2)
%>  out1 = arg2;
%>  out2 = arg1;
%>end
%>@endverbatim
%>


% Add HiveTemp directories to the search path list of Octave
if ~exist('HT_VAR_LIB_PATH')
  error('Library path is not defined <HT_VAR_LIB_PATH>');
endif

if ~is_absolute_filename(HT_VAR_LIB_PATH)
  HT_VAR_LIB_PATH = make_absolute_filename(HT_VAR_LIB_PATH);
endif

addpath(strcat(HT_VAR_LIB_PATH, '/Files'));
addpath(strcat(HT_VAR_LIB_PATH, '/Math'));
addpath(strcat(HT_VAR_LIB_PATH, '/Tools'));

% Load miscellaneous package for functions "clip", "normc"
pkg load miscellaneous;
pkg load geometry;      % boundingBox
pkg load linear-algebra; % rotv

% Import constant symbols that could be useful to design a thermal model
HT_ImportConstants();

% ======================= Global variables definitions =========================
% Instead of specifying explicit argument with function calls, it is slightly
% more efficient to use global variables for very specific test

%> @var HT_VAR_EPSILON_POS
%> @brief It defines the smallest meaningful variation of position
%>      With <1E-6> it means no object can be smaller than 1 micrometer. This
%>      variable allows perfect position comparison to be done otherwise, there
%>      are always issues with numeric rounding.
global HT_VAR_EPSILON_POS;
HT_VAR_EPSILON_POS = 1E-6;

%> @var It is similar to @ref HT_VAR_EPSILON_POS but for normalized
%>     position. It defines the smallest meaningful variation of normalized
%>     position used mostly by faces objects. The ratio of @ref HT_VAR_EPSILON_POS and
%>     @ref HT_VAR_EPSILON_U defines the biggest object than can be handle by the
%>     library. Example: with 1E-6 and 1E-8 (default values), the biggest object
%>     would be 100m before numeric rounding error may cause problems
global HT_VAR_EPSILON_U;
HT_VAR_EPSILON_U = 1E-8;

%> @var It defines if input parameters must be checked during function startup
%>     This variable is not yet taken into account in all functions but should
%>     be in the future
global HT_VAR_CHECK_INPUT;
HT_VAR_CHECK_INPUT = true;

%> @var It defines if unused or unknown input parameters must be checked during
%>     function startup. If true, functions will raise an error if an unknown
%>     parameter name is specified.
global HT_VAR_STRICT_INPUT_PARAMETERS;
HT_VAR_STRICT_INPUT_PARAMETERS = true;



