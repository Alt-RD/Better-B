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
% Build a 1-dimensionnal multilayer heat conduction problem along the X axis
% Thickness of first and last nodes are halved.
%
% Input arguments:
% .name = [string] nom du model
% .params = [struct] data structure containing following fields
%         .length = [1xNlayer] length of each layer along the x direction
%         .material  = struct array/cell array of {material}
%         .rc     = [1x(Nlayer-1)] thermal contact resistance
%         .dim    = 2x1 wall dimension along y and z directions
%         .axis   = [3x3]=[X,Y,Z] axis of local coordinate system. It defines
%                   the wall orientation
%         .globalPosition = [3x1] position of the wall (x,y,z)=(0,0,0) in
%                   global coordinate system.
%         .n      = [1xNlayer] number of nodes in each layer
%                      (Nlayer >= 2) (extreme nodes are located on boundaries)
%         .T0     = [1x1] ou [1xnNode] initial temperature defined globally or
%                   for each layer
%         .mergeFaces = [logical] if true laterals faces are merged into on faces
%                       when a multilayer volume is defined
%         .gridType = 'ff'/'hh'/'hf'/'fh' defines the type of boundary nodes
%                     "full" or "half" nodes
%         .boundaryNames : Replace the name of first and last node by XX.boundaryname(1) and XX.boundaryname(2)
% .options = [struct] data structure containing following options:
%         .verbose
%         .nodeNameModel = [string] model name used to specify a name for each
%                          node. By default: 'n%d'. <%d> will be replaced by
%                          node index.
%                          If two arguments '%d' are defined, the first will be
%                          replaced by layer index, and second by node index
%                          relative to current layer.
%
% Output arguments:
% .M = [data structure of type 'model_conduction_1d']. See "HT_BuildModel" for a list of existing fields
% .faces = array of data stucture of type 'face'.
%          Local face axis are defined as follow. It shares the same axis with the model
%          and are ordered in such a way the third vector is pointing outside the
%          volume.
%          Face orientation is ordered based on local axis: -x,+x,-y,+y,-z,+z
% .nodes = [cell array of string] list of node names defined in this model
function [M faces nodes volumes] = HT_Model_Conduction1D(name, varargin)
  HT_ImportConstants();

  lTicId = tic();
  lModelTypeName = 'model_conduction_1D';

  assert(~isempty(strtrim(name)), 'Invalid model name. It is empty.');

  lParams = struct();
  options = struct();

  if isstruct(varargin{1})        % If first parameter is a struct, parameters are loaded from it
    lParams = varargin{1};

    if (numel(varargin) > 1) && isstruct(varargin{end})
      options = varargin{end};
      prop = varargin(2:2:(end-1));
      values = varargin(3:2:(end-1));
    else
      prop = varargin(2:2:end);
      values = varargin(3:2:end);
    endif
  else
    if isstruct(varargin{end})
      options = varargin{end};
      prop = varargin(1:2:(end-1));
      values = varargin(2:2:(end-1));
    else
      options = struct();
      prop = varargin(1:2:end);
      values = varargin(2:2:end);
    endif
  endif

  % Set default parameters for structure <options>
  options = CheckField(options, 'verbose',        true,   { @(v) islogical(v) } );
  options = CheckField(options, 'verboselevel',   true,   { @(v) isscalar(v) && (v >= 0) && (round(v) == v) } );
  options = CheckField(options, 'skipFieldCheck', false,    @(v) islogical(v));
  if options.verbose, disp(sprintf('Building models <%s> type <%s>', name, lModelTypeName)); endif;

  % Add to the <lParams> structure all parameters specified
  for i=1:numel(prop)
    lParams = setfield(lParams, prop{i}, values{i});
  endfor

  % Check parameter structure fields. It is desirable to check that no silly errors
  % are made on the parameter name
  if ~options.skipFieldCheck
    lValidParamList = { 'mergeFaces', 'length', 'base', 'n', 'dim', 'globalPosition', ...
                        'gridType', 'material', 'rc', 'T0', 'axis', 'nodeNameModel', 'boundaryNames' };
    lParamList = fieldnames(lParams);
    lValidParam = cellfun(@(v) any(strcmp(v, lValidParamList)), lParamList);
    assert(all(lValidParam), sprintf('Invalid parameter name detected <%s>', strjoin(lParamList(~lValidParam), ',')));
    clear lParamList lValidParam lValidParamList;
  endif

  % Set default parameter values
  lParams = CheckField(lParams, 'material',       [],                 {"exist", @(v) isstruct(v) && iscolumn(v) });
  lParams = CheckField(lParams, 'n',              [],                 {"exist", @(v) (isfloat(v) || iscell(v)) && iscolumn(v) });
  lParams = CheckField(lParams, 'u',              [],                 @(v) (isfloat(v) || iscell(v)) && iscolumn(v));
  lParams = CheckField(lParams, 'T0',             NaN,                @(v) (isfloat(v) || iscell(v)) && iscolumn(v));
  lParams = CheckField(lParams, 'axis',           [],                 @(v) isfloat(v) && all(size(v) == [3,3]));
  lParams = CheckField(lParams, 'globalPosition', [],                 @(v) isfloat(v) && (numel(v) == 3));
  lParams = CheckField(lParams, 'length',         [],                 {"exist", @(v) (isscalar(v) || iscell(v)) && iscolumn(v) });
  lParams = CheckField(lParams, 'dim',            [],                 @(v) isfloat(v) || (numel(v) == 2) );
  lParams = CheckField(lParams, 'base',           [],                 @(v) HT_CheckType(v, 'face') );
  lParams = CheckField(lParams, 'nodeNameModel',  strcat(strrep(name, ' ', '_'), '.n%d'),       @(v) ischar(v) && any(numel(strfind(options.nodeNameModel, '%d')) == [ 1 2 ] ));
  lParams = CheckField(lParams, 'boundaryNames',  { {'',''} },        @(v) iscell(v));
  lParams = CheckField(lParams, 'gridType',       'ff',               @(v) ischar(v) || (iscell(v) && iscolumn(v)));
  lParams = CheckField(lParams, 'mergeFaces',     false,              @(v) islogical(v));

  % Build a structure array of each layer
  % If an argument has an invalid size, Octave will raise an error here
  lLayerData = struct(
    'material', arrayfun(@(v) v, lParams.material, 'UniformOutput', false), ...
    'n', lParams.n, ...
    'u', lParams.u, ...     % Position of node boundaries
    'upos', [], ...         % Position of nodes
    'usize', [], ...        % Size of nodes
    'T0', lParams.T0, ...
    'length', lParams.length, ...
    'boundaryNames', lParams.boundaryNames,...
    'gridType', lParams.gridType,...
    'V', [], ...
    'C', [], ...
    'RVec', [], ...
    'G', [],...
    'nodes', []);

  % The number of layer is equal to the number of element in the structure array
  nLayer = numel(lLayerData);

  % Check parameters for each layer
  for i=1:numel(lLayerData)
    layer = lLayerData(i);
    layer = CheckField(layer, 'material',           [],         @(v) HT_CheckType(v, 'material') );
    layer = CheckField(layer, 'n',                  [],         @(v) isfloat(v) && isscalar(v) && (round(v)==v) && (v > 0));
    layer = CheckField(layer, 'u',                  [],         @(v) isempty(v) || (isnumeric(v) && (numel(v) == layer.n) && all(v > 0) && (u(1) == 0) && (u(end) == 1)));
    layer = CheckField(layer, 'T0',                 [],         @(v) isfloat(v) && any(numel(v) == [1 layer.n]));
    layer = CheckField(layer, 'length',             [],         @(v) isscalar(v) && (v > 0));
    layer = CheckField(layer, 'boundaryNames',      {'',''},    @(v) iscell(v) && (numel(v) == 2) && all( cellfun(@(x) ischar(x), v )) );
    layer = CheckField(layer, 'gridType',           'ff',       @(v) ischar(v) && any(strcmpi(v, {'ff', 'fh', 'hf', 'hh'})));

    % <u> contains the nodes boundaries
    if isempty(layer.u) % If not specified by the user, regularly spaced nodes is used
      ldx = 1 / (layer.n - 0.5*sum(layer.gridType == 'h'));
      layer.usize = repmat(ldx, layer.n, 1);
      if (layer.gridType(1) == 'h'), layer.usize(1) *= 0.5; endif
      if (layer.gridType(2) == 'h'), layer.usize(end) *= 0.5; endif
      layer.u = [0; cumsum(layer.usize)];
      clear ldx;
    else
      layer.usize = layer.u(2:end) - layer.u(1:(end-1));
    endif

    layer.upos = 0.5*(layer.u(1:(end-1)) + layer.u(2:end));
    if (layer.gridType(1) == 'h'), layer.upos(1) = 0; endif
    if (layer.gridType(2) == 'h'), layer.upos(end) = 1; endif

    if isscalar(layer.T0), layer.T0 = repmat(layer.T0, 1, layer.n); endif

    lLayerData(i) = layer;
  endfor

  lParams = CheckField(lParams, 'rc',             zeros(nLayer-1,1),  @(v) isnumeric(v) && any(numel(v) == [1, nLayer-1]) && all(v >= 0));

  % Build the function that construct node names depending on
  % content of <nodeNameMode>
  if numel(strfind(lParams.nodeNameModel, '%d')) == 1
    NodeNameFunc = @(layer, node, totalnode) sprintf(lParams.nodeNameModel, totalnode);
  else
    NodeNameFunc = @(layer, node, totalnode) sprintf(lParams.nodeNameModel, layer, node);
  endif

  % The field "base" is used to specify a base to the new model.
  % It is generally a face on which the 1D conduction model is built
  if ~isempty(lParams.base)
    assert(HT_CheckType(lParams.base, 'face'), 'Invalid structure provided for field <base>');
    assert(isempty(lParams.dim) && isempty(lParams.globalPosition), 'Field <dim> and <globalPosition> can not be specified when a base is provided');

    % The axis is extracted from the face
    lFace = lParams.base;
    lParams.axis = cross(lFace.norm, lFace.axis(:,U)); % A direct coordinate system is built using cross product between Normal vector and U
    lParams.axis = [lFace.norm lFace.axis(:,U) lParams.axis];

    % If a coordinate system is provided, the coordinate system of the specified face
    % is changed accordingly.
    lFace = HT_Face_AdjustAxis(lFace, lParams.axis);

    lParams.dim = lFace.size;
    lParams.globalPosition = lFace.globalPosition;
    clear lFace;
  else
    % Set default value of <axis> and <globalPosition>
    if isempty(lParams.dim),            lParams.dim = ones(2,1); endif
    if isempty(lParams.axis),           lParams.axis = eye(3); endif
    if isempty(lParams.globalPosition), lParams.globalPosition = zeros(3,1); endif
  endif

  nTotal = sum(arrayfun(@(v) v.n, lLayerData)); % Total number of nodes
  lArea = prod(lParams.dim); % Surface du mur dans le plan xy

  % Create matrix for each layer
  lNodeCount = 0;
  for i=1:nLayer
    layer = lLayerData(i);
    layer.V = (lArea * layer.length) * layer.usize';
    layer.C = (lArea * layer.length * HT_Material_GetRhoC(layer.material)) * layer.usize';
    layer.C = layer.V * HT_Material_GetRhoC(layer.material);
    layer.nodes = cell(1, layer.n);

    for k=1:layer.n
      layer.nodes{k} = NodeNameFunc(i, k, lNodeCount+k);
    endfor

    % lRVec contains the thermal resistance between each nodes and boundaries
    layer.RVec = [layer.upos - layer.u(1:(end-1)); 0] + [0; layer.u(2:end) - layer.upos];
    layer.RVec *= (layer.length / HT_Material_GetLambda(layer.material, X));
    % Create the conductance matrix
    % G(i,i) = -S/R(i,i-1)-S/R(i,i+1)
    GVec = lArea./layer.RVec;
    GVec(1) = 0;      % Extreme nodes are not inserted now in the matrix G
    GVec(end) = 0;
    layer.G = spdiags([GVec(2:end), -(GVec(1:(end-1)) + GVec(2:end)), GVec(1:(end-1))], [-1 0 1], layer.n, layer.n);

    if ~isempty(layer.boundaryNames{1})
      layer.nodes{1}    = sprintf(layer.boundaryNames{1}, strrep(name, ' ', '_'));
    endif
    if ~isempty(layer.boundaryNames{2}) && (layer.n > 1)
      layer.nodes{end}  = sprintf(layer.boundaryNames{2}, strrep(name, ' ', '_'));
    endif

    lNodeCount += layer.n;
    lLayerData(i) = layer;
  endfor
  clear GVec layer;

  % Build the model
  M = HT_Model_Init(lModelTypeName, name, lParams);
##  M.nodes = HT_CellUnwrap({lLayerData.nodes}); [lLayerData](:);   % Extract and merge all nodes from all layers
  M.nodes = [lLayerData.nodes](:);   % Extract and merge all nodes from all layers
  M.G = blkdiag(lLayerData.G);    % Extract and merge all conductance matrix
  M.C = [lLayerData.C](:);
  M.T0 = [lLayerData.T0](:);
  M.globalPosition = lParams.globalPosition;
  M.axis = lParams.axis;

  if isargout(1)
    % Add contact resistance and thermal connection between layer
    lNodeCount = 0;
    for i=1:(nLayer-1)
      layer = lLayerData(i);
      lG = lArea / (layer.RVec(end) + lParams.rc(i) + lLayerData(i+1).RVec(1));
      ind = lNodeCount+layer.n;
      M.G(ind, ind) -= lG;
      M.G(ind+1, ind+1) -= lG;
      M.G(ind, ind+1) += lG;
      M.G(ind+1, ind) += lG;
      lNodeCount += layer.n;
    endfor

    assert(all(abs(sum(M.G,2)) < 1E-12), 'Invalid conductance matrix'); % Vérifie que chaque ligne a une somme nulle
  endif

  if isargout(2) % Build face array
    ly = lParams.dim(1);
    lz = lParams.dim(2);

    if lParams.mergeFaces
      lLambdaVec = arrayfun(@(v) HT_Material_GetLambda(v.material), lLayerData); % Thermal conductivity of each layer
      lLengthVec = arrayfun(@(v) v.length, lLayerData);
      lNVec = arrayfun(@(v) v.n, lLayerData);
      lx = sum(lLengthVec);

      faces = HT_Face_Init({ 'Xm', 'Xp', 'Ym', 'Yp', 'Zm', 'Zp'}, ...
                              'size', { [lz ly], [ly lz], [lx lz], [lz lx], [ly lx], [lx ly]}, ...
                              'globalPosition', { lParams.globalPosition, ...
                                                  lParams.globalPosition + lx * lParams.axis(:,X), ...
                                                  lParams.globalPosition, ...
                                                  lParams.globalPosition + ly * lParams.axis(:,Y), ...
                                                  lParams.globalPosition, ...
                                                  lParams.globalPosition + lz * lParams.axis(:,Z)}, ...
                              'n', { [1 1]; [1 1]; [nTotal 1]; [1 nTotal]; [1 nTotal]; [nTotal 1] },...
                              'axis', {[lParams.axis(:,3) lParams.axis(:,2)], ... zy
                                       [lParams.axis(:,2) lParams.axis(:,3)], ... yz
                                       [lParams.axis(:,1) lParams.axis(:,3)], ... xz
                                       [lParams.axis(:,3) lParams.axis(:,1)], ... zx
                                       [lParams.axis(:,2) lParams.axis(:,1)], ... yx
                                       [lParams.axis(:,1) lParams.axis(:,2)] }, ... xy;
                              'norm', { -lParams.axis(:,X), ...  %Face XM oriented to outside
                                        lParams.axis(:,X), ...
                                        -lParams.axis(:,Y), ...
                                        lParams.axis(:,Y), ...
                                        -lParams.axis(:,Z), ...
                                        lParams.axis(:,Z) }, ...
                              'r', {  lLayerData(1).RVec(1), ...
                                      lLayerData(end).RVec(end), ...
                                      ly/2./lLambdaVec,... .*dx_vec*lz),... 10/01/2023: resistance is a density not absolute value
                                      ly/2./lLambdaVec,... .*dx_vec*lz),...
                                      lz/2./lLambdaVec,... .*dx_vec*ly),...
                                      lz/2./lLambdaVec},... .*dx_vec*ly)}, ...
                              'material', lParams.material, ...
                              'model', name);

      % Pour chaque face, les axes sont les mêmes que le repère principal 3D
      % L'ordre est tel que les axes forment toujours un repère direct.
      % Ex: Sur la face YP, le repère 2D est (u,v)=(z3D,x3D)
      % Ex: Sur la face ZP, le repère 2D est (u,v)=(x3D,y3D)

      faces(FaceXM).nodes = M.nodes(1);
      faces(FaceXM).pos = [0.5  0.5]; % Axis y x
      faces(FaceXM).dims = [0 1 0 1]; ...

      faces(FaceXP).nodes = M.nodes(nTotal);
      faces(FaceXP).pos = [0.5  0.5]; % Axis x y
      faces(FaceXP).dims = [0 1 0 1]; ...

      % Calcul de la position et de la taille de chaque noeud
      lLengthVec_cumsum = cumsum([0; lLengthVec(1:(end-1))]);
      lLengthVec_rep = repelem(lLengthVec_cumsum, lNVec);

      % Position of each nodes in real space
      P = cell2mat(arrayfun(@(v) v.upos * v.length, lLayerData, 'UniformOutput', false));
      P += lLengthVec_rep;
      P /= lx;

      P1 = cell2mat(arrayfun(@(v) v.u(1:(end-1)) * v.length, lLayerData, 'UniformOutput', false));
      P1 += lLengthVec_rep;
      P1 /= lx;
      P2 = cell2mat(arrayfun(@(v) v.u(2:end) * v.length, lLayerData, 'UniformOutput', false));
      P2 += lLengthVec_rep;
      P2 /= lx;

      faces(FaceYM).nodes = M.nodes(1:nTotal);
      faces(FaceYM).pos = [P repmat(0.5, nTotal,1)]; % u,v = Axis x z
      faces(FaceYM).dims = [P1 , P2 , ...
                        zeros(nTotal, 1) , ones(nTotal, 1)];

      faces(FaceYP).nodes = M.nodes(1:nTotal);
      faces(FaceYP).pos = [repmat(0.5, nTotal,1) P]; % u,v = Axis z x
      faces(FaceYP).dims = [zeros(nTotal, 1)  , ones(nTotal, 1) ...
                        P1            , P2];

      faces(FaceZM).nodes = M.nodes(1:nTotal);
      faces(FaceZM).pos = [repmat(0.5, nTotal,1) P]; % u,v = y,x
      faces(FaceZM).dims = [zeros(nTotal, 1)  ,  ones(nTotal, 1) , ...
                        P1            ,  P2];

      faces(FaceZP).nodes = M.nodes(1:nTotal);
      faces(FaceZP).pos = [P repmat(0.5, nTotal,1)]; % u,v = x,y
      faces(FaceZP).dims = [P1            , P2, ...
                        zeros(nTotal, 1)  , ones(nTotal, 1)];
    else
##      lLengthVec = arrayfun(@(v) v.length, lLayerData);
      lLengthVec = [lLayerData.length](:);
      lLengthVec_cumsum = cumsum([0; lLengthVec(1:(end-1))]);
      lNVec = arrayfun(@(v) v.n, lLayerData);
      lNCumSumVec = cumsum([0; lNVec]);
##      lMaterialIndexList = repelem((1:nLayer)', lNVec);
      faces = [];

      for i=1:nLayer
        layer = lLayerData(i);
        lx = layer.length;
        lLambda = HT_Material_GetLambda(layer.material);

        f = HT_Face_Init({ 'Xm'; 'Xp'; 'Ym'; 'Yp'; 'Zm'; 'Zp'}, ...
                              'size', { [lz ly]; [ly lz]; [lx lz]; [lz lx]; [ly lx]; [lx ly]}, ...
                              'globalPosition', { lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X); ...
                                                  lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X) + lx * lParams.axis(:,X); ...
                                                  lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X); ...
                                                  lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X) + ly * lParams.axis(:,Y); ...
                                                  lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X); ...
                                                  lParams.globalPosition + lLengthVec_cumsum(i) * lParams.axis(:,X) + lz * lParams.axis(:,Z)}, ...
                              'n', { [1 1]; [1 1]; [layer.n 1]; [1 layer.n]; [1 layer.n]; [layer.n 1] },...
                              'axis', {[lParams.axis(:,3) lParams.axis(:,2)]; ... zy
                                       [lParams.axis(:,2) lParams.axis(:,3)]; ... yz
                                       [lParams.axis(:,1) lParams.axis(:,3)]; ... xz
                                       [lParams.axis(:,3) lParams.axis(:,1)]; ... zx
                                       [lParams.axis(:,2) lParams.axis(:,1)]; ... yx
                                       [lParams.axis(:,1) lParams.axis(:,2)] }, ... xy;
                              'norm', { -lParams.axis(:,X); ...  %Face XM oriented to outside
                                        lParams.axis(:,X); ...
                                        -lParams.axis(:,Y); ...
                                        lParams.axis(:,Y); ...
                                        -lParams.axis(:,Z); ...
                                        lParams.axis(:,Z) }, ...
                              'r', {  lLayerData(i).RVec(1); ...
                                      lLayerData(i).RVec(end); ...
                                      ly/2/lLambda(Y);...
                                      ly/2/lLambda(Y);...
                                      lz/2/lLambda(Z);...
                                      lz/2/lLambda(Z)},...
                              'material', {layer.material; ... Only one material per face here
                                           layer.material;
                                           layer.material; ...
                                           layer.material; ...
                                           layer.material; ...
                                           layer.material}, ...
                              'nodes', {M.nodes(lNCumSumVec(i)+1); ...
                                        M.nodes(lNCumSumVec(i)+layer.n);
                                        M.nodes(lNCumSumVec(i)+(1:layer.n)); ...
                                        M.nodes(lNCumSumVec(i)+(1:layer.n)); ...
                                        M.nodes(lNCumSumVec(i)+(1:layer.n)); ...
                                        M.nodes(lNCumSumVec(i)+(1:layer.n))}, ...
                              'materialIndex', {1;
                                                1;
                                                ones(layer.n, 1);
                                                ones(layer.n, 1);
                                                ones(layer.n, 1);
                                                ones(layer.n, 1)}, ...
                              'model', name);

        % Pour chaque face, les axes sont les mêmes que le repère principal 3D
        % L'ordre est tel que les axes forment toujours un repère direct.
        % Ex: Sur la face YP, le repère 2D est (u,v)=(z3D,x3D)
        % Ex: Sur la face ZP, le repère 2D est (u,v)=(x3D,y3D)
        f(FaceXM).pos = [0.5  0.5]; % Axis y x
        f(FaceXM).dims = [0 1 0 1]; ...

        f(FaceXP).pos = [0.5  0.5]; % Axis x y
        f(FaceXP).dims = [0 1 0 1]; ...

        % Position of each nodes in real space
        P = layer.upos;
        P1 = layer.u(1:(end-1));
        P2 = layer.u(2:end);

        f(FaceYM).pos = [P repmat(0.5, layer.n,1)]; % u,v = Axis x z
        f(FaceYM).dims = [P1 , P2 , ...
                          zeros(layer.n, 1) , ones(layer.n, 1)];

        f(FaceYP).pos = [repmat(0.5, layer.n,1) P]; % u,v = Axis z x
        f(FaceYP).dims = [zeros(layer.n, 1)  , ones(layer.n, 1) ...
                          P1            , P2];

        f(FaceZM).pos = [repmat(0.5, layer.n,1) P]; % u,v = y,x
        f(FaceZM).dims = [zeros(layer.n, 1)  ,  ones(layer.n, 1) , ...
                          P1            ,  P2];

        f(FaceZP).pos = [P repmat(0.5, layer.n,1)]; % u,v = x,y
        f(FaceZP).dims = [P1            , P2, ...
                          zeros(layer.n, 1)  , ones(layer.n, 1)];

        faces = [faces ; f];
      endfor
    endif
  endif

  if isargout(4)
    volumes = HT_Volume_Init('name', arrayfun(@(layerIndex) sprintf('%s.layer%d', name, layerIndex), 1:nLayer, 'UniformOutput', false), ...
                             'shapeType', 'cube', ...
                             'nodes', arrayfun(@(layer) layer.nodes, lLayerData, 'UniformOutput', false), ...
                             'nodesVolume', arrayfun(@(layer) layer.V, lLayerData, 'UniformOutput', false));
  endif

  if nargout > 2
    nodes = M.nodes;
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
  else  % 12/01/2023: put the function test in the "else" condition to allow invalid data to be specified as default
    for i=1:numel(chkFunc)
      if strcmp(class(chkFunc{i}), 'function_handle')
        assert(chkFunc{i}(getfield(P, field)), sprintf('Check failed: field <%s>', field));
      endif
    endfor
  endif

endfunction

