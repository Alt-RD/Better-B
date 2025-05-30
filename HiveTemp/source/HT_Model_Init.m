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
% This function returns an empty model structure
%
% Input arguments:
% <type> = [string] or {cell of string} specifying the model type
% <name> = [string] the name of the model
% <params> = [struct] data stucture ('paramname', paramvalue, ...) 
%            containing parameters
%
function M = HT_Model_Init(type, name, params)
  assert(nargin >= 1, 'model type must be defined');
  
  if nargin < 2, name = ''; endif;
  if nargin < 3, params = struct(); endif;
  
  if (iscell(type)) % Add the type 'model'
    if size(type,2) > 0
      type = unique([type, 'model']);
    else
      type = unique([type; 'model']);
    endif
  elseif (ischar(type))
    type = {'model', type};
  else
    error('Argument <type> is not valid');
  endif

  M = struct( 'name', name, ...
              '__type__', [],...
              'globalPosition', [], ...
              'G', sparse(0,0), ...
              'C', [], ...
              'Gdyn', [], ...         % Time varying coefficient
              'nodes', [], ...
              'T0', [], ...           % Initial temperature
              'rad', [], ...          % Radiation models
              'params', params, ...
              'axis', [], ...
              'submodel', [], ...
              'timer', [], ...
              'counter', []);
              
  M.__type__ = type;  % Type may be a cell, since a structure may have severals types
  M.rad = {}; % Can not initialize with an empty cell since struct interpet this as a struct array
  M.timer = struct('buildingTime', 0, 'solvingTime', 0);
  M.counter = struct('warnings', 0, 'errors', 0);
endfunction
