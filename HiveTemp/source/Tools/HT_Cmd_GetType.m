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
function infos = HT_Cmd_GetType(cmd, options)
  if nargin < 2, options = struct(); endif;

  if strcmpi(cmd.type, 'temperature')
    infos = INT_GetTemperatureType(cmd, options);
  elseif strcmpi(cmd.type, 'flux')
    infos = INT_GetFluxType(cmd, options);
  else
    error('Invalid command type');
  endif
endfunction

function infos = INT_GetTemperatureType(cmd, options)
  infos = struct( 'type', 0, ...
                  'userData', false, ...
                  'typeString', '');

  if isnumeric(cmd.data)
    infos.type = INT_GetTemperatureDataType(cmd.data, numel(cmd.nodes), cmd.name);
    infos.typeString = sprintf('Cmd <%s> is <%s>', cmd.name, HT_Cmd_GetTypeName(infos.type));

  elseif iscell(cmd.data)
    % If data is a cell array, it means face were specified in the field <nodes>
    assert(iscell(cmd.nodes), 'Invalid type of <data>');
    assert(numel(cmd.data) == numel(cmd.nodes), 'Invalid size of <data>');

    if numel(cmd.data) == 1
      infos.type = INT_GetTemperatureDataType(cmd.data{1}, numel(cmd.nodes{1}.nodes), cmd.name);
      infos.typeString = sprintf('Cmd <%s> is <%s> (single face)', cmd.name, HT_Cmd_GetTypeName(infos.type));
    else
      infos.type = zeros(numel(cmd.data), 1);
      for k=1:numel(cmd.data)
        infos.type(k) = INT_GetTemperatureDataType(cmd.data{k}, numel(cmd.nodes{k}.nodes), cmd.name);
      endfor
      infos.typeString = sprintf('Cmd <%s> is multi type <%s> (multi faces)', cmd.name, sprintf('%d,', infos.type));
    endif

  else
    error('Not implemented');
##    assert(is_function_handle(cmd.data), 'Command must be matrix or function handle');
##    lArgs = Int_GetFunctionArguments(cmd.data);
##    % DO guess the command type, the last argument <userData> if present is removed
##    if any(strcmpi(lArgs{end}, {'user', 'userData', 'data'}))
##      lArgs(end) = [];
##      infos.userData = true;
##    endif
##    % If cmd.nodes were a cell array of face, the <face> argument must be the last one
##    if cmd.sys.faceProvided
##      assert(strcmpi(lArgs{end}, 'face'), ...
##      'Argument <face> of <.data> function must be present, and before the optional <user> argument');
##      lArgs(end) = [];
##    endif
##
##    assert(numel(lArgs) == 1);
##    if strcmpi(lArgs{1}, 't')
##      infos.type = 15; % Function with time in argument   @(t)
##      infos.typeString = sprintf('  -> Temperature type cmd <%s> is <Variable (function)(t) and uniform temperature>', cmd.name);
##    elseif numel(lArgs) == 1
##      infos.type = 14; % Function with index in argument  @(k) or an other name
##      infos.typeString = sprintf('  -> Temperature type cmd <%s> is <Variable (function)(index) and uniform temperature>', cmd.name);
##    else
##      error('Invalid argument for function <.data>');
##    endif
  endif
endfunction


function T = INT_GetTemperatureDataType(data, nNodes, name)
  HT_ImportConstants();

  if columns(data) == 2         % Interpolation array
    T = HT_BTYPE_T_INTERPOLATION_ARRAY; % Interpolated
  elseif isscalar(data)         % Constant and uniform temperature
    T = HT_BTYPE_T_CONSTANT_UNIFORM;
  elseif iscolumn(data) && (numel(data) == nNodes)        % Constant temperatuer specified for each nodes
    T = HT_BTYPE_T_CONSTANT_NON_UNIFORM;
  elseif isrow(data)            % Variable and uniform temperature
    T = HT_BTYPE_T_VARIABLE_UNIFORM;
  elseif (rows(data) == nNodes) && (columns(data) > 2)
    T = HT_BTYPE_T_VARIABLE_NON_UNIFORM;
  else
    error('Cmd_GetType::Invalid command data');
  endif

endfunction





function infos = INT_GetFluxType(cmd, options)
  HT_ImportConstants();

  infos = struct( 'type', 0, ...
                  'userData', false, ...
                  'typeString', '');

  lUnit = 'undefined';

  if iscell(cmd.area)
    lAreaValidFlag = cellfun(@(a) ~isempty(a), cmd.area);

    if all(~lAreaValidFlag)
      lUnit = 'absolute';
    elseif all(lAreaValidFlag)
      lUnit = 'density';
    else
      lUnit = 'mixte';
    endif
  else
    error('Not implemented');
  endif

  if isnumeric(cmd.data)
    infos.type = INT_GetFluxDataType(cmd.data, numel(cmd.nodes), cmd.name);
    infos.typeString = sprintf('Cmd <%s> is <%s> [unit=%s]', cmd.name, HT_Cmd_GetTypeName(infos.type), lUnit);

  elseif iscell(cmd.data)
    % If data is a cell array, it means face were specified in the field <nodes>
    assert(iscell(cmd.nodes), 'Invalid type of <data>');
    assert(numel(cmd.data) == numel(cmd.nodes), 'Invalid size of <data>');
    % Check that every data has the same type

    if numel(cmd.data) == 1
      infos.type = INT_GetFluxDataType(cmd.data{1}, numel(cmd.nodes{1}.nodes), cmd.name);
      infos.typeString = sprintf('Cmd <%s> is <%s> (single face) [unit=%s]', cmd.name, HT_Cmd_GetTypeName(infos.type), lUnit);
    else
      infos.type = zeros(numel(cmd.data), 1);
      for k=1:numel(cmd.data)
        infos.type(k) = INT_GetFluxDataType(cmd.data{k}, numel(cmd.nodes{k}.nodes), cmd.name);
      endfor
##      infos.typeString = sprintf('Cmd <%s> is multi type <%s> (multi faces =%d)', cmd.name, sprintf('%d,', infos.type), numel(cmd.data));

      lFaceFlag = cellfun(@(f) HT_CheckType(f, 'face'), cmd.nodes);
      lFaceNodeCount = sum(cellfun(@(f) numel(f.nodes), cmd.nodes(lFaceFlag)));

      infos.typeString = sprintf('Cmd <%s> is multi type <%s> (%d faces, %d nodes, %d total nodes) [unit=%s]', ...
                                            cmd.name, ...
                                            HT_Cmd_GetTypeName(unique(infos.type)), ...
                                            sum(lFaceFlag), ...
                                            sum(~lFaceFlag), ...
                                            sum(~lFaceFlag) + lFaceNodeCount, ...
                                            lUnit);
    endif

  else
    error('Not implemented');
##    assert(is_function_handle(cmd.data), 'Command must be matrix or function handle');
##    lArgs = Int_GetFunctionArguments(cmd.data);
##    % DO guess the command type, the last argument <userData> if present is removed
##    if any(strcmpi(lArgs{end}, {'user', 'userData', 'data'}))
##      lArgs(end) = [];
##      infos.userData = true;
##    endif
##    % If cmd.nodes were a cell array of face, the <face> argument must be the last one
##    if cmd.sys.faceProvided
##      assert(strcmpi(lArgs{end}, 'face'), ...
##      'Argument <face> of <.data> function must be present, and before the optional <user> argument');
##      lArgs(end) = [];
##    endif
##
##    assert(any(numel(lArgs) == [1 3]), sprintf('Field <%s.data> is not valid. Invalid argument count.', cmd.name));
##
##    if (numel(lArgs) == 1) && strcmpi(lArgs{1}, 't')
##      infos.type = 27; % Function with time in argument   @(t)
##    elseif (numel(lArgs) == 1)
##      infos.type = 26; % Function with index in argument  @(k) or an other name
##    elseif (numel(lArgs) == 3)
##      infos.type = 28; % Function @(t,k,index)
##      if ~strcmpi(lArgs{1}, 't')
##        warning(sprintf('First argument of function <%s.data> should be called "t"', cmd.name));
##      endif
##    else
##      error('Invalid argument for function <.data>.');
##    endif
  endif

endfunction



function T = INT_GetFluxDataType(data, nNodes, name)
  HT_ImportConstants();

  if columns(data) == (nNodes+1)                  % Interpolation array
    T = HT_BTYPE_FLUX_INTERPOLATION_ARRAY; % Interpolated

  elseif isscalar(data)                           % Constant flux and uniform
    T = HT_BTYPE_FLUX_CONSTANT_UNIFORM;

  elseif isrow(data)                              % Variable flux and uniform
    T = HT_BTYPE_FLUX_VARIABLE_UNIFORM;

  elseif iscolumn(data)                              % Constant flux and non uniform
    assert(numel(cmd.data) == nNodes);
    T = HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM;

  elseif (rows(data) == nNodes) && (columns(data) > 1) % Variable flux and non uniform
    T = HT_BTYPE_FLUX_VARIABLE_NON_UNIFORM;

  else
    error(sprintf('Field <%s.data> is not valid', name));
  endif

endfunction


% Check the type of a function handle
% time varying or non linear time varying
% Returns cell(nfunction, 1) containing a array of str
function P = Int_GetFunctionArguments(func)
  assert(is_function_handle(func));

  funStr = functions(func).function;
  ind1 = strfind(funStr, '(');
  ind2 = strfind(funStr, ')');
  P = strsplit(funStr((ind1+1):(ind2-1)), {',', ' '});
endfunction
