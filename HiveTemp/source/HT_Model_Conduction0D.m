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
% Build a 0-dimensionnal multilayer heat conduction problem
%
% Input arguments:
% .name = [string] nom du model
% .params = [struct] data structure containing following fields
%         .rhoC   = [1xNlayer] volumetric thermal capacity of each layer
%         .lambda = [1xNlayer] thermal conductivity of each layer
%         .dim    = 3x1 volume dimensions
%         .axis   = [3x3]=[X,Y,Z] axis of local coordinate system. It defines
%                   the volume orientation 
%         .globalPosition = [3x1] position of the volume (x,y,z)=(0,0,0) in
%                   global coordinate system.
%         .T0     = [1x1] initial temperature
% .options = [struct] data structure containing following options:
%         .verbose
%         .nodeName = [string] model node name to replace default name
%
% Output arguments:
% .M = [data structure of type 'model_conduction_0d']. See "HT_Model_Init" for a list of existing fields
% .faces = array of data stucture of type 'face'.
%          Local face axis are defined as follow. It shares the same axis with the model
%          and are ordered in such a way the third vector is pointing outside the
%          volume.
function [M faces] = HT_Model_Conduction0D(name, params, options)
  lTicId = tic();
  lModelTypeName = 'model_conduction_0D';
  
  if options.verbose, disp(sprintf('Building models <%s> type <%s>', name, lModelTypeName)); tic; endif; % Start timer (tic)
    
  assert(nargin >= 2, 'There must be at least 2 input arguments');
  
  if nargin < 3, options = struct(); endif;
  assert(numel(options) == 1, 'Invalid <option> struct');
  assert(numel(params) == 1, 'Invalid <params> struct');
  
  % Set default parameters for structure <options>
  options = CheckField(options, 'verbose', true, {@(v) islogical(v)});
  options = CheckField(options, 'nodeName', '');
  
  % Set default parameters for structure <params>
  params = CheckField(params, 'dim', [], {"exist", @(v) any(numel(v)==[1 3]) && all(v > 0)});
  
  if isfield(params, 'material')
    assert(HT_CheckType(params.material, "material"), 'Invalid material structure type');
    assert(~isfield(params, 'rhoC') && ~isfield(params, 'lambda'), 'Material parameter conflict');
    params.rhoC = HT_Material_GetRhoC(params.material);
    params.lambda = HT_Material_GetLambda(params.material);
  endif
  
  params = CheckField(params, 'rhoC', [],   {"exist", @(v) all(v > 0) && numel(v) == 1});
  params = CheckField(params, 'lambda', [], {"exist", @(v) all(v > 0) && any(numel(v) == [1 3])});
  params = CheckField(params, 'dim', [], {"exist", @(v) any(numel(v)== [1 3]) && all(v > 0)});
  params = CheckField(params, 'T0', NaN, @(v) (numel(v) == 1));

  if ~isfield(params, 'axis') && options.verbose
    disp(sprintf('Warning: axis are not defined for model <%s>', name));
  endif
  
  if numel(params.lambda) == 1, params.lambda = ones(3,1) * params.lambda; endif;
  
  params = CheckField(params, 'axis', eye(3),   {@(v) all(size(v) == [3 3]) && all(vecnorm(v, Inf, 1) == 1.0)});
  
  if ~isfield(params, 'globalPosition') && options.verbose
    disp(sprintf('Warning: global position is not defined for model <%s>', name));
  endif

  params = CheckField(params, 'globalPosition', zeros(3,1),   {@(v) all(size(v) == [3 1])});
  
  M = HT_Model_Init(lModelTypeName, name, params);
    
  % Default value for nodeNameModel
  if isempty(options.nodeName), options.nodeName = strcat(strrep(name, ' ', '_'), '.n0'); endif;
  
  % Count the number of occurence of %d  
  idx = strfind(options.nodeName, '%');
  assert(isempty(idx), 'nodeName can not contains special char %%');
  
  % Construction des matrices
  M.nodes = cell(1, 1);
  M.nodes(1) = options.nodeName;
  M.G = sparse(1,1);
  M.C = prod(params.dim) * params.rhoC; % Initialisation du vecteur des capacités
  M.T0 = params.T0;
  M.globalPosition = params.globalPosition;
  M.axis = params.axis;
  
  if nargout > 1 % Build face array
    lx = params.dim(1);
    ly = params.dim(2);
    lz = params.dim(3);
    
    X=1; Y=2; Z=3;
    
    fclip = @(v,l,b) min(max(v,l), b);

    faces = HT_Face_Init({ 'Xm', 'Xp', 'Ym', 'Yp', 'Zm', 'Zp'}, ...
                            'size', { [lz ly], [ly lz], [lx lz], [lz lx], [ly lx], [lx ly]}, ...
                            'globalPosition', { M.globalPosition, ...
                                                M.globalPosition + lx * params.axis(:,1), ...
                                                M.globalPosition, ...
                                                M.globalPosition + ly * params.axis(:,2), ...
                                                M.globalPosition, ...
                                                M.globalPosition + lz * params.axis(:,3)}, ...
                            'axis', {[M.axis(:,3) M.axis(:,2)], ... zy
                                     [M.axis(:,2) M.axis(:,3)], ... yz
                                     [M.axis(:,1) M.axis(:,3)], ... xz
                                     [M.axis(:,3) M.axis(:,1)], ... zx
                                     [M.axis(:,2) M.axis(:,1)], ... yx
                                     [M.axis(:,1) M.axis(:,2)] }, ... xy;
                            'norm', { -params.axis(:,1), ...  %Face XM oriented to outside
                                      params.axis(:,1), ...
                                      -params.axis(:,2), ...
                                      params.axis(:,2), ...
                                      -params.axis(:,3), ...
                                      params.axis(:,3) }, ...
                            'r', {  lx/2./(params.lambda(X)), ...
                                    lx/2./(params.lambda(X)), ...
                                    ly/2./(params.lambda(Y)),...
                                    ly/2./(params.lambda(Y)),...
                                    lz/2./(params.lambda(Z)),...
                                    lz/2./(params.lambda(Z))}, ...
                            'model', name);
    
    XM=1; XP=2; YM=3; YP=4; ZM=5; ZP=6;
    % Pour chaque face, les axes sont les mêmes que le repère principal 3D
    % L'ordre est tel que les axes forment toujours un repère direct.
    % Ex: Sur la face YP, le repère 2D est (u,v)=(z3D,x3D)
    % Ex: Sur la face ZP, le repère 2D est (u,v)=(x3D,y3D)
    
    faces(XM).nodes = M.nodes(1);
    faces(XM).pos = [0.5 0.5]; % z, y
    faces(XM).dims = [0 1 0 1];

    faces(XP).nodes = M.nodes(1);
    faces(XP).pos = [0.5 0.5]; % z, y
    faces(XP).dims = [0 1 0 1];
    
    faces(YM).nodes = M.nodes(1);
    faces(YM).pos = [0.5 0.5]; % z, y
    faces(YM).dims = [0 1 0 1];
    
    faces(YP).nodes = M.nodes(1);
    faces(YP).pos = [0.5 0.5]; % z, y
    faces(YP).dims = [0 1 0 1];
    
    faces(ZM).nodes = M.nodes(1);
    faces(ZM).pos = [0.5  0.5]; % Axis y x
    faces(ZM).dims = [0 1 0 1]; ...

    faces(ZP).nodes = M.nodes(1);
    faces(ZP).pos = [0.5  0.5]; % Axis x y
    faces(ZP).dims = [0 1 0 1]; ...
 
  endif

  M.timer.buildingTime += toc(lTicId);
  
  if options.verbose, disp(sprintf('Building models <%s> done in %.2f ms', name, 1000*M.timer.buildingTime)); endif; % Start timer (tic)

endfunction


% Vérification des paramètres d'entrée
function P = CheckField(P, field, value, chkFunc)
  assert(nargin >= 3);
  
  if nargin < 4, chkFunc = {}; endif;
  
  % Convert to cell array
  if ~iscell(chkFunc), chkFunc = {chkFunc}; endif;
  
  lChkExist = false;
    
  for i=1:numel(chkFunc)
    if strcmpi(chkFunc{i}, 'exist')
      lChkExist = true;
      break;
    endif
  endfor
  
  if ~isfield(P, field)
    assert(~lChkExist, sprintf('Field <%s> does not exist', field));
    
    P = setfield(P, field, value);
  elseif ~iscell(getfield(P, field))
    P = setfield(P, field, cast(getfield(P, field), class(value)));
  endif
  
  for i=1:numel(chkFunc)
    if strcmp(class(chkFunc{i}), 'function_handle')
      assert(chkFunc{i}(getfield(P, field)), sprintf('Check failed: field <%s>', field));
    endif
  endfor
  
endfunction
  
