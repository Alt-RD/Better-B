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
%
% Construit un modèle de conduction 3D multilayer
% Arguments:
% .name = nom du model
% .nodeNameModel = 'n%d_%d_%d' chaine avec l'indice %d des noeuds
%                  Par défaut: <nom_du_model>.n%d_%d_%d
%
% Retourne une structure model
% .name = nom du model
% .G = NxN sparse matrice du modèle
% .C
% .nodes = {1xdim N} liste des noeuds
% .Gdyn
% .T0
% .type = "model_conduction_1d"
% .data = params
%
% La structure <params> contient les champs
% .dim = [2x1] longueur directions x,y
% .length = [1xNlayer] longueur de chaque couche dans la direction z
% .rhoC = [1xNlayer] capacité thermique de chaque couche
% .lambda = [3xNlayer] conductivités thermiques de chaque couche
% .Rc =  [1x(Nlayer-1)] résistance de contact
% .n = [NX x NY x NZ1...xNZ_Nlayer] discrétisation de chaque couche (nombre de noeuds)
%                 (les noeuds extrêmes sont sur les bords)
% .xgrid = [1 x NX] spécifie la taille du maillage selon x  (normalized)
%          [NX x 1] spécifie la position des lignes du maillage  (normalized) (0 is implicit, not specified)
%          { pos(1), n(1), pos(2), n(2), ... pos(i), n(i) }  (normalized)
% .xgridtype = 'hh'/'hf'/'fh'/'ff' -> node type at XmXp position: 'h' for half node and 'f' for full node
% .ygrid = [1 x NY] spécifie la taille du maillage selon y   (normalized)
%          [NY x 1] spécifie la position des lignes du maillage  (normalized)
%          { pos(1), n(1), pos(2), n(2), ... pos(i), n(i) }
% .ygridtype = 'hh'/'hf'/'fh'/'ff'
% .zgrid = [1 x NY] spécifie la taille du maillage selon z   (normalized)
%          [NY x 1] spécifie la position des lignes du maillage  (normalized)
%          { pos(1), n(1), pos(2), n(2), ... pos(i), n(i) } (normalized)
%          Coordinate pos(i) may be several time 1 indicating end of current layer
%          Otherwise, it is normalized by the full material length.
% .zgridtype = 'hh'/'hf'/'fh'/'ff' -> Only for extreme node (not between layers)
% .T0 = température initiale [1x1] ou [1xnNode]
% .axis = dim3x3 = [x y z] composantes des vecteurs x,y,z dans le repère global (défaut: Identity)
% .mergeFaces = [true/false]
%
% La structure <options> contient les champs suivant:
% .verbose
% .verboselevel
% .nodeNameModel : model name for nodes
% .skipFieldCheck : [true/false] prevent the parameter name check to be done
%
% .faces : retourne les 6 faces dans une structure face
% 'Xm', 'Xp', 'Ym', 'Yp', 'Zm', 'Zp'
%       .nodes
%       .pos
%       .dims
%       .name
%       .type
%       .size
function [M, faces] = HT_Model_Conduction3D(name, params, options)
  HT_ImportConstants();

  lTicId = tic();
  lModelTypeName = 'model_conduction_3D';

  assert(nargin >= 2, 'There must be at least 2 input arguments');

  if nargin < 3, options = struct(); endif;
  assert(numel(options) == 1, 'Invalid <option> struct');
  assert(numel(params) == 1, 'Invalid <params> struct');

  options = HT_CheckField(options, 'verbose', true,            @(v) islogical(v));
  options = HT_CheckField(options, 'verboselevel',   true,     { @(v) isscalar(v) && (v >= 0) && (round(v) == v) } );
  options = HT_CheckField(options, 'nodeNameModel', '');
  options = HT_CheckField(options, 'skipFieldCheck', false,    @(v) islogical(v));

  if options.verbose, disp(sprintf('Building models <%s> type <%s>', name, lModelTypeName)); endif;

  % Check parameter structure fields. It is desirable to check that no silly errors
  % are made on the parameter name
  if ~options.skipFieldCheck
    lValidParamList = { 'mergeFaces', 'length', 'base', 'n', 'dim', 'globalPosition', ...
                       'xGrid', 'xGridType', 'yGrid', 'yGridType', 'zGrid', 'zGridType', 'material', 'rhoC', 'lambda', 'Rc', 'T0', 'axis' };
    lParamList = fieldnames(params);
    lValidParam = cellfun(@(v) any(strcmp(v, lValidParamList)), lParamList);
    assert(all(lValidParam), sprintf('Invalid parameter name detected <%s>', strjoin(lParamList(~lValidParam), ',')));
    clear lParamList lValidParam lValidParamList;
  endif

  params = HT_CheckField(params, 'mergeFaces',       true, @(v) islogical(v));
  params = HT_CheckField(params, 'length',           [], {"exist", @(v) all(v > 0)});
  params = HT_CheckField(params, 'base',             [], @(v) HT_CheckType(params.base, 'face'));
  params = HT_CheckField(params, 'n',                NA(3,1));
  params = HT_CheckField(params, 'axis',             [], @(v) isempty(v) || all(size(v)==[3 3]));
  params = HT_CheckField(params, 'dim',              NA(2,1));
  params = HT_CheckField(params, 'globalPosition',   NA(3,1));

  % The field "base" is used to specify a base to the new model.
  % It is generally a face on which the 3D conduction model is built
  if ~isempty(params.base)
    %assert(~isfield(params, 'dim') || isempty(params.dim), 'Field <dim> can not be specified when a base is provided');
    %assert(~isfield(params, 'globalPosition') || isempty(params.globalPosition), 'Field <globalPosition> can not be specified when a base is provided');
    assert(isempty(params.axis) || xor(all(isna(params.axis)), all(~isna(params.axis))), 'Field <axis> is not valid');

    % The axis is extracted from the face
    lFace = params.base;
    % A direct coordinate system is built using cross product between Normal vector and U
    % If the user already specified an axis, it is kept
    if isempty(params.axis) || all(isna(params.axis))
      params.axis = [lFace.axis(:,U), cross(lFace.norm, lFace.axis(:,U)), lFace.norm];
    endif

    % Adapt the coordinate system of the face to the new local coordinate system
    % The local axis is extracted for x,y direction of local system x,y,z specified by the user
    % z is the normal and has nothing to do here
    lFace = HT_Face_AdjustAxis(lFace, params.axis(:,[X Y]));
    % If not specified by the user, dimension vector is set
    if isempty(params.dim), params.dim = NA(2,1); endif
    params.dim = postpad(params.dim, 2, NA);
    params.dim(isna(params.dim)) = lFace.size(isna(params.dim));

    % If not specified by the user, dimension vector is set
    if isempty(params.globalPosition), params.globalPosition = NA(3,1); endif
    params.globalPosition = postpad(params.globalPosition, 3, NA);
    params.globalPosition(isna(params.globalPosition)) = lFace.globalPosition(isna(params.globalPosition));

    % Overwrite xGrid and xGridType only if not specified by the user
    % If the user specified node count <.n>, nothing is loaded from the <base> face
    if (~isfield(params, 'xGrid') || isempty(params.xGrid)) && (~isfield(params, 'n') || isempty(params.n) || isna(params.n(X)))
      % yGrid and zGrid is also retrieved from the base face
      if isna(params.n(X))
        params.n(X) = lFace.n(U);
      endif

      params.xGrid = unique(lFace.dims(:,2));
      % If nodes of the base face are on the boundary, it means we must use 'half' node mode

      if ~isfield(params, 'xGridType') || isempty(params.xGridType)
        params.xGridType = 'ff';
        if abs(min(lFace.pos(:,U)))   < 1E-10, params.xGridType(1) = 'h'; endif
        if abs(max(lFace.pos(:,U))-1) < 1E-10, params.xGridType(2) = 'h'; endif
      endif
    endif

    % Overwrite yGrid and yGridType only if not specified by the user
    if (~isfield(params, 'yGrid') || isempty(params.yGrid)) && (~isfield(params, 'n') || isempty(params.n) || isna(params.n(Y)))
      if isna(params.n(Y))
        params.n(Y) = lFace.n(V);
      endif

      params.yGrid = unique(lFace.dims(:,4));

      % If nodes of the base face are on the boundary, it means we must use 'half' node mode
      if ~isfield(params, 'yGridType') || isempty(params.yGridType)
        params.yGridType = 'ff';
        if abs(min(lFace.pos(:,V)))   < 1E-10, params.yGridType(1) = 'h'; endif
        if abs(max(lFace.pos(:,V))-1) < 1E-10, params.yGridType(2) = 'h'; endif
      endif
    endif

    clear lFace;
  endif

  Nlayer = numel(params.length);
  lz = sum(params.length);
  params = HT_CheckField(params, 'n', [], {"exist", @(v) isfloat(v) && all((v > 0) & (round(v) == v) | isna(v) ) && any(numel(v) == [1 Nlayer]+2)});

  % It allows a single face or a structure array of faces to be specified in xGrid field.
  % (not only cell array)
  % It avoid the use of double {{ }} that would be necessary
  % params.xGrid = {{ lMod_Underroof_faces{FaceZP} }}...
  if isfield(params, 'xGrid') && isstruct(params.xGrid)
    lFaceList = params.xGrid;
    params.xGrid = cell(numel(lFaceList), 1);
    for i=1:numel(lFaceList)
      params.xGrid{i} = lFaceList(i);
    endfor
  endif

  if isfield(params, 'yGrid') && isstruct(params.yGrid)
    lFaceList = params.yGrid;
    params.yGrid = cell(numel(lFaceList), 1);
    for i=1:numel(lFaceList)
      params.yGrid{i} = lFaceList(i);
    endfor
  endif

  if isfield(params, 'zGrid') && isstruct(params.zGrid)
    lFaceList = params.zGrid;
    params.zGrid = cell(numel(lFaceList), 1);
    for i=1:numel(lFaceList)
      params.zGrid{i} = lFaceList(i);
    endfor
  endif

  params = HT_CheckField(params, 'dim',        [], {'exist', @(v) isfloat(v) && all(v > 0) && (numel(v) == 2) && all(~isna(v))});
  params = HT_CheckField(params, 'xGrid',      [], {@(v) isempty(v) || (numel(v) == params.n(X)) || (iscell(v) && mod(numel(v),2) == 0) || (iscell(v) && all(cellfun(@(x) HT_CheckType(x, 'face'), v))) });
  params = HT_CheckField(params, 'yGrid',      [], {@(v) isempty(v) || (numel(v) == params.n(Y)) || (iscell(v) && mod(numel(v),2) == 0) || (iscell(v) && all(cellfun(@(x) HT_CheckType(x, 'face'), v))) });
  params = HT_CheckField(params, 'xGridType',  'hh', {@(v) any(strcmpi(v, {'hh', 'hf', 'fh', 'ff'}))});
  params = HT_CheckField(params, 'yGridType',  'hh', {@(v) any(strcmpi(v, {'hh', 'hf', 'fh', 'ff'}))});

  params = HT_CheckField(params, 'zGrid',      [], {@(v) isempty(v) || (numel(v) == params.n(Z)) || ...
                                                 (iscell(v) && (mod(size(v, 2),2) == 0) && any(size(v,1) == [1 Nlayer]) )});
  params = HT_CheckField(params, 'zGridType',  'hh', {@(v) any(strcmpi(v, {'hh', 'hf', 'fh', 'ff'}))});

  params = HT_CheckField(params, 'material',   [], @(v) HT_CheckType(params.material, 'material') );
  if ~isempty(params.material)
    params.rhoC = HT_Material_GetRhoC(params.material);
    params.lambda = HT_Material_GetLambda(params.material);
    if isscalar(params.lambda), params.lambda = repmat(params.lambda, 3, Nlayer); endif
  endif

  params = HT_CheckField(params, 'rhoC',       [], {'exist', @(v) isfloat(v) && all(v > 0) && any(numel(v) == [1 Nlayer])});
  params = HT_CheckField(params, 'lambda',     [], {'exist', @(v) isfloat(v) && all(v > 0) && (size(v,1) == 3) && (size(v, 2) == Nlayer)});
  params = HT_CheckField(params, 'Rc', zeros(Nlayer-1,1), @(v) isfloat(v) && all(v >= 0) && (any(numel(v) == [1 Nlayer-1])));
  params = HT_CheckField(params, 'T0', NaN, @(v) (numel(v) == 1) || (numel(v) == nNodes));

  params = HT_CheckField(params, 'axis', eye(3),   {@(v) all(size(v) == [3 3]) && all(vecnorm(v, Inf, 1) == 1.0)});
  assert(norm(cross(params.axis(:,1), params.axis(:,2)) - params.axis(:,3)) < 1E-15, 'Invalid <axis>');
  params = HT_CheckField(params, 'globalPosition', zeros(3,1),   {@(v) (numel(v) == 3) && all(~isna(v))});

  assert(isna(params.n(X)) || (params.n(X) > 1) || strcmpi(params.xGridType, 'ff'), sprintf('xGridType <%s> is not compatible with 1 node only', params.xGridType));
  assert(isna(params.n(Y)) || (params.n(Y) > 1) || strcmpi(params.yGridType, 'ff'), sprintf('yGridType <%s> is not compatible with 1 node only', params.yGridType));
  assert(isna(params.n(Z)) || (params.n(Z) > 1) || strcmpi(params.zGridType(1), 'f'), sprintf('zGridType <%s> is not compatible with 1 node only', params.zGridType));
  assert(isna(params.n(end)) || (params.n(end) > 1) || strcmpi(params.zGridType(2), 'f'), sprintf('zGridType <%s> is not compatible with 1 node only', params.zGridType));

  % Position, size, and boundaries of nodes
  [xPos, xSize, xDim, params.n(1)]      = Int_InitMeshSize(params.xGrid, params.n(1), params.dim(1), params.xGridType, params.axis(:,X), params.globalPosition);
  [yPos, ySize, yDim, params.n(2)]      = Int_InitMeshSize(params.yGrid, params.n(2), params.dim(2), params.yGridType, params.axis(:,Y), params.globalPosition);
  [zPos, zSize, zDim, params.n(3:end)]  = Int_InitMeshSize(params.zGrid, params.n(3:end), params.length, params.zGridType, params.axis(:,Z), params.globalPosition);
  dxVec = (xPos(2:end) - xPos(1:(end-1))) * params.dim(1);
  dyVec = (yPos(2:end) - yPos(1:(end-1))) * params.dim(2);
  dzVec = (zPos(2:end) - zPos(1:(end-1))) * lz;

  assert(~any(isna(params.n)), 'Invalid mesh. Some node count are not defined');

  M = HT_Model_Init(lModelTypeName, name, params);

  % Default value for nodeNameModel
  if isempty(options.nodeNameModel), options.nodeNameModel = strcat(strrep(name, ' ', '_'), '.n%d'); endif;

  nx = params.n(1);
  ny = params.n(2);
  nzVec = params.n(3:(2+Nlayer));
  nz = sum(nzVec);
  nxy = params.n(1)*params.n(2);
  nNodes = nxy * nz;

  % Count the number of occurence of %d
  idx = numel(strfind(options.nodeNameModel, '%d'));
  assert(any(idx == [1 2 3 4]), 'nodeNameModel must contains at least 3 occurences of <%d>');

  if idx == 1,
    NodeNameFunc = @(layer, knode, i,j,k) sprintf(options.nodeNameModel, (k-1)*nxy + (j-1)*nx + i);
  elseif idx == 2,
    NodeNameFunc = @(layer, knode, i,j,k) sprintf(options.nodeNameModel, layer, (knode-1)*nxy + (j-1)*nx + i);
  elseif idx == 3,
    NodeNameFunc = @(layer, knode, i,j,k) sprintf(options.nodeNameModel, i,j,k);
  elseif idx == 4,
    NodeNameFunc = @(layer, knode, i,j,k) sprintf(options.nodeNameModel, layer, i,j,knode);
  endif
  clear idx;

  if numel(params.rhoC) == 1, params.rhoC = repmat(params.rhoC, nNodes, 1); endif;
  if numel(params.lambda) == 1, params.lambda = repmat(params.lambda, nNodes, 1); endif;
  if numel(params.T0) == 1, params.T0 = repmat(params.T0, nNodes, 1); endif;
  if numel(params.Rc) == 1, params.Rc = repmat(params.Rc, Nlayer-1, 1); endif;

  % Construction des matrices
  Gmat = sparse(nNodes,nNodes);
  nodes = cell(nNodes, 1);

  % Vecteur des capacités et matrice des conductances
  Cvec = repmat(params.dim(1)*params.dim(2)*lz, 1);
  Cvec .*= repmat(xSize, ny*nz, 1);
  Cvec .*= repmat(repelem(ySize, nx, 1), nz, 1);
  Cvec .*= repelem(zSize, nx*ny, 1);

  xSize *= params.dim(1);
  ySize *= params.dim(2);
  zSize *= lz;

  iNode = 0;
  for m=1:Nlayer
    nzlayer = params.n(2+m);

    ind = sum(params.n(2+(1:(m-1)))) * nxy;
    Cvec((ind+1):(ind + nzlayer*nxy)) *= params.rhoC(m);

    klayer = sum(params.n(3:(2+m-1))); % Indice k du début de la couche

    lambdax = params.lambda(1, m);
    lambday = params.lambda(2, m);
    lambdaz = params.lambda(3, m);

    for k=1:nzlayer
      for j=1:ny
        for i=1:nx
          iNode++;
          ktotal = klayer + k;
          nodes{iNode} = NodeNameFunc(m, k, i,j, ktotal);

          Gx = lambdax* zSize(ktotal)*ySize(j);
          Gy = lambday* xSize(i)*zSize(ktotal);
          Gz = lambdaz* xSize(i)*ySize(j);

          % Conductances
          if i < nx
            Gmat(iNode, iNode+1) += Gx/dxVec(i);
            Gmat(iNode, iNode) -= Gx/dxVec(i);
          endif
          if i > 1
            Gmat(iNode, iNode-1) += Gx/dxVec(i-1);
            Gmat(iNode, iNode) -= Gx/dxVec(i-1);
          endif

          if j < ny
            Gmat(iNode, iNode+nx) += Gy/dyVec(j);
            Gmat(iNode, iNode) -= Gy/dyVec(j);
          endif
          if j > 1
            Gmat(iNode, iNode-nx) += Gy/dyVec(j-1);
            Gmat(iNode, iNode) -= Gy/dyVec(j-1);
          endif

          % Next node in z direction
          if k < nzlayer
            Gmat(iNode, iNode+nxy) += Gz/dzVec(klayer+k);
            Gmat(iNode, iNode) -= Gz/dzVec(klayer+k);
          elseif m < Nlayer % Connect with next layer
            Gzp = params.lambda(3, m+1) * Sz / zSize(ktotal+1);
            Gzc = lambdaz * Sz / zSize(ktotal);
            lGz = 1/(0.5/Gzc + 0.5/Gzp);
            Gmat(iNode, iNode+nxy) += lGz;
            Gmat(iNode, iNode) -= lGz;
            clear lGz;
          endif
          if k > 1
            Gmat(iNode, iNode-nxy) += Gz/dzVec(klayer+k-1);
            Gmat(iNode, iNode) -= Gz/dzVec(klayer+k-1);
          elseif m > 1 % Connect with previous layer
            Gzp = params.lambda(3, m-1) * Sz / zSize(ktotal-1);
            Gzc = lambdaz * Sz / zSize(ktotal);
            lGz = 1/(0.5/Gzc + 0.5/Gzp);
            Gmat(iNode, iNode-nxy) += lGz;
            Gmat(iNode, iNode) -= lGz;
            clear lGz;
          endif
        endfor
      endfor
    endfor

  endfor

  M.nodes = nodes;
  M.G = Gmat;
  M.C = Cvec; % Initialisation du vecteur des capacités
  M.T0 = params.T0;
  M.axis = params.axis;
  M.globalPosition = params.globalPosition(:);

  assert(all(abs(sum(Gmat,2) ./ max(Gmat, [], 2)) < 1E-14)); % Vérifie que chaque ligne a une somme nulle
  assert(Gmat == Gmat', 'Conductance matrix is not diagonal');

  clear nodes Gmat iNode Cvec lambdax lambday lambdaz Sx Sy Sz;

  if nargout > 1
    if params.mergeFaces
      lx = params.dim(1);
      ly = params.dim(2);
      lambdaxVec = [];
      lambdayVec = [];
      for m=1:Nlayer
        lambdaxVec = [lambdaxVec; repmat(params.lambda(1,m), nzVec(m), 1)];
        lambdayVec = [lambdayVec; repmat(params.lambda(2,m), nzVec(m), 1)];
      endfor

      faces = HT_Face_Init({ 'Xm'; 'Xp'; 'Ym'; 'Yp'; 'Zm'; 'Zp'}, ...
                              'size', { [lz ly]; [ly lz]; [lx lz]; [lz lx]; [ly lx]; [lx ly]}, ...
                              'globalPosition', { M.globalPosition; ...
                                                  M.globalPosition + params.dim(1) * params.axis(:,1); ...
                                                  M.globalPosition; ...
                                                  M.globalPosition + params.dim(2) * params.axis(:,2); ...
                                                  M.globalPosition; ...
                                                  M.globalPosition + lz * params.axis(:,3)}, ...
                              'n', { [nz ny]; [ny nz]; [nx nz]; [nz nx]; [ny nx]; [nx ny] },...
                              'axis', {[M.axis(:,3) M.axis(:,2)]; ... zy
                                       [M.axis(:,2) M.axis(:,3)]; ... yz
                                       [M.axis(:,1) M.axis(:,3)]; ... xz
                                       [M.axis(:,3) M.axis(:,1)]; ... zx
                                       [M.axis(:,2) M.axis(:,1)]; ... yx
                                       [M.axis(:,1) M.axis(:,2)] }, ... xy;
                              'norm', { -params.axis(:,1); ...  %Face XM oriented to outside
                                        params.axis(:,1); ...
                                        -params.axis(:,2); ...
                                        params.axis(:,2); ...
                                        -params.axis(:,3); ...
                                        params.axis(:,3) }, ...
                               'material', params.material, ...
                               'model', name...
                                       );

      % Pour chaque face, les axes sont les mêmes que le repère principal 3D
      % L'ordre est tel que les axes forment toujours un repère direct.
      % Ex: Sur la face YP, le repère 2D est (u,v)=(z3D,x3D)
      % Ex: Sur la face ZP, le repère 2D est (u,v)=(x3D,y3D)

      % Make sure nodes are following axis order (y first, then x)
      faces(FaceZM).nodes = M.nodes( reshape(1:nxy, nx, ny)'(:) );
      faces(FaceZM).pos =  [repmat(yPos, nx, 1)   repelem(xPos, ny, 1)]; % Axis y x
      faces(FaceZM).dims = [repmat(yDim, nx, 1)   repelem(xDim, ny, 1)];
      if params.zGridType(1) != 'h', faces(FaceZM).r = 0.5*zSize(1)/params.lambda(3,1); endif;

      faces(FaceZP).nodes = M.nodes( (nz-1)*nxy + (1:nxy) );
      faces(FaceZP).pos =  [repmat(xPos, ny, 1)   repelem(yPos, nx, 1)]; % Axis x y
      faces(FaceZP).dims = [repmat(xDim, ny, 1)   repelem(yDim, nx, 1)];
      if params.zGridType(2) != 'h', faces(FaceZP).r = 0.5*zSize(end)/params.lambda(3,end); endif;

      faces(FaceXM).nodes = M.nodes( 1 + reshape((0:(ny*nz-1))*nx, ny, nz)'(:) );
      faces(FaceXM).pos =  [repmat(zPos, ny, 1)    repelem(yPos, nz, 1)]; % z, y
      faces(FaceXM).dims = [repmat(zDim, ny, 1)    repelem(yDim, nz, 1)];
      if params.xGridType(1) != 'h'
        faces(FaceXM).r = 0.5*xSize(1)./ repmat(lambdaxVec, ny, 1); % .* repmat(zSize, ny, 1) .* repelem(ySize, nz, 1) );
      endif;

      faces(FaceXP).nodes = M.nodes( reshape((1:(ny*nz))*nx, ny, nz)(:) );
      faces(FaceXP).pos =  [repmat(yPos, nz, 1)     repelem(zPos, ny, 1)];  % Axis y z
      faces(FaceXP).dims = [repmat(yDim, nz, 1)     repelem(zDim, ny, 1)];
      if params.xGridType(2) != 'h'
        faces(FaceXP).r = 0.5*xSize(end)./ repelem(lambdaxVec, ny, 1); % .* repelem(zSize, ny, 1) .* repmat(ySize, nz, 1) );
      endif;

      faces(FaceYM).nodes = M.nodes( repmat((1:nx)', 1, nz)(:) + repmat((0:(nz-1))*nxy, nx, 1)(:) );
      faces(FaceYM).pos =  [repmat(xPos, nz, 1)  repelem(zPos, nx, 1)]; % Axis x z
      faces(FaceYM).dims = [repmat(xDim, nz, 1)  repelem(zDim, nx, 1)];
      if params.yGridType(1) != 'h'
        faces(FaceYM).r = 0.5*ySize(1)./ repelem(lambdayVec, nx, 1); %  .* repmat(xSize, nz, 1) .* repelem(zSize, nx, 1) );
      endif;

      faces(FaceYP).nodes = M.nodes( repmat((nxy-nx+1):nxy, nz, 1)(:) + repmat((0:(nz-1))'*nxy, nx, 1) );
      faces(FaceYP).pos =  [repmat(zPos, nx, 1)  repelem(xPos, nz, 1)];    % Axis z x
      faces(FaceYP).dims = [repmat(zDim, nx, 1)  repelem(xDim, nz, 1)];
      if params.yGridType(2) != 'h'
        faces(FaceYP).r = 0.5*ySize(end)./ repmat(lambdayVec, nx, 1); %  .* repmat(zSize, nx, 1) .* repelem(xSize, nz, 1) );
      endif;
    else
      error('Not implemented');
    endif
  endif

  M.timer.buildingTime += toc(lTicId);

  if options.verbose, disp(sprintf('Building models <%s> done in %.2f ms', name, 1000*M.timer.buildingTime)); endif; % Start timer (tic)

endfunction


% nx [Nlayer x 1] number of nodes per layer
% Lx [Nlayer x 1] length of each layer
% xtype = 'xx' with x=h/f for node type at boundaries
% vDir [dim 3x1] direction in global axis
% vPos [dim 3x1] global position
% Returns:
% xPos = node position (normalized)
% xSize = node size (normalized)
% xDim = [Nx1] [xstart xend; ...] with N=sum(nx)
% nx = [Nx1] the number of nodes for each layer
function [xPos xSize xDim, nx] = Int_InitMeshSize(xGrid, nx, Lx, xtype, vDir, vPos)
  assert(numel(nx) == numel(Lx), 'Invalid size of node count vector. Must be equal to the number of layers.');
  Nlayer = numel(nx);
  Ltx = sum(Lx);
  Ltx_cum = cumsum([0; Lx]);
  ntx = sum(nx);

  xPos = [];
  xSize = [];

  if isempty(xGrid)
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    fdx = @(m) Lx(m) / Ltx / (nx(m) - 0.5*((m==1) && (xtype(1)=='h')) - 0.5*((m==Nlayer) && (xtype(2)=='h')));

    for i=1:Nlayer
      ldx = fdx(i);
      lxSize = repmat(ldx, nx(i), 1);
      if i == 1 && (xtype(1)=='h'), lxSize(1) *= 0.5; endif; % Half node ?
      if i == Nlayer && (xtype(2)=='h'), lxSize(end) *= 0.5; endif; % Half node ?

      lxPos = sum(Lx(1:(i-1))/Ltx) + ((1:nx(i))'-0.5)*ldx - (0.5*ldx * ((xtype(1)=='h') && (i==1)));

      xSize = [xSize ; lxSize];
      xPos = [xPos ; lxPos];
    endfor

    assert(abs(xPos(end)-1+0.5*ldx*(xtype(1)!='h')) < 1E-12);
    assert(abs(sum(xSize)-1) < 1E-12);
    clear ldx lxSize lxPos;
  elseif iscell(xGrid) && all(cellfun(@(v) isnumeric(v), xGrid))
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    % Converts coordinates in xGrid vector
    if size(xGrid, 2) > 1, xGrid = xGrid'; endif;

    xGridpos = cell2mat(xGrid);
    xGridn = xGridpos(2:2:end);
    xGridpos = xGridpos(1:2:end);
    assert(all(xGridn > 0), 'Invalid node count in xyzgrid vector');
    assert(all(xGridpos > 0), 'Invalid node position in xyzgrid vector');

    lOneInd = find(abs(xGridpos-1) < 1E-14);
    if numel(lOneInd) == Nlayer % Ex: { 0.0, 4, 0.5, 4, 1.0, 10, 1.0, 8 } % 1.0 indicates next layer
      for i=1:numel(lOneInd)
        if i == 1
          xGridpos(1:lOneInd(i)) *= Lx(i) / Ltx;
          assert(sum(xGridn(1:lOneInd(i))) == nx(i), 'Wrong node count');
        else
          xGridpos((lOneInd(i-1)+1):lOneInd(i)) = xGridpos(lOneInd(i-1)) + ...
                                                  xGridpos((lOneInd(i-1)+1):lOneInd(i)) * Lx(i)/Ltx;
          assert(sum(xGridn((lOneInd(i-1)+1):lOneInd(i))) == nx(i), 'Wrong node count');
        endif

      endfor

      clear lOneInd;
    else
      assert(numel(lOneInd) == 1, 'Pos 1.0 should appear once (or nlayer)');
      xGridpos = xGrid(1:2:end);
      xGridn = xGrid(2:2:end);
    endif

    assert(issorted(xGridpos));

    % xGrid contient les positions relatives et les nombres de noeuds
    fdx = @(L, m, mtotal) L / (xGridn(m) - 0.5*((m==1) && (xtype(1)=='h')) - 0.5*((m==mtotal) && (xtype(2)=='h')));

    % Cas particulier pour le premier intervalle
    nPart = numel(xGridpos);
    ldx = fdx(xGridpos(1), 1, nPart);
    lxSize = repmat(ldx, xGridn(1), 1);
    if xtype(1)=='h', lxSize(1) *= 0.5; endif; % Half node ?
    if ((nPart==1) && (xtype(2)=='h')), lxSize(end) *= 0.5; endif; % Half node ?
    xSize = lxSize;
    xPos = ((1:xGridn(1))'-0.5)*ldx - (lxSize(1) * (xtype(1)=='h'));

    for i=2:nPart
      ldx = fdx(xGridpos(i)-xGridpos(i-1), i, nPart);
      lxSize = repmat(ldx, xGridn(i), 1);
      if i == nPart && (xtype(2)=='h'), lxSize(end) *= 0.5; endif; % Half node ?

      lxPos = xGridpos(i-1) + ((1:xGridn(i))'-0.5)*ldx;

      xSize = [xSize ; lxSize];
      xPos = [xPos ; lxPos];
    endfor

    assert(abs(xPos(end)-1+0.5*ldx*(xtype(1)!='h')) < 1E-12);
    assert(abs(sum(xSize) -1) < 1E-12);
    clear ldx lxSize lxPos nPart;
  elseif iscell(xGrid) && all(cellfun(@(v) HT_CheckType(v, 'face'), xGrid))
% nx [Nlayer x 1] number of nodes per layer
% Lx [Nlayer x 1] length of each layer
% xtype = 'xx' with x=h/f for node type at boundaries
% Returns:
% xPos = node position (normalized)
% xSize = node size (normalized)
% xDim = [Nx1] [xstart xend; ...] with N=sum(nx)
    lFaceList = xGrid;
    lFaceCount = numel(lFaceList);
    lFaceGrid = [];
    for i=1:lFaceCount
      F = lFaceList{i};
      F = HT_Face_AdjustAxis(F, vDir);
      lGridVertex = HT_Face_GetAbsoluteData(F, 'xGrid');
      % Compute the dot product between each vertex (vertex - model.globalPosition) and grid direction
      lGridPos = vDir' * (lGridVertex - vPos);
      % And convert the row vector to column vector
      lFaceGrid = [lFaceGrid; lGridPos'];
    endfor

    lFaceGrid /= Ltx; % normalized

    % Remove invalid values (if face is higher than current model)
    lFaceGrid((lFaceGrid < 0) | (lFaceGrid > 1)) = [];
    % Now add the layers
    lFaceGrid = [lFaceGrid; 0; Lx/Ltx];
    lFaceGrid = unique(lFaceGrid); % Unique and sorted

    % The final number of nodes is either the number specified by the user
    % or the number infered by the list of faces.
    lNodeCount = max(ntx, numel(lFaceGrid)-1);

    if Nlayer == 1
      % If the number of nodes was not specified
      % <xN> is set and returned to the main function
      if isna(ntx), ntx = lNodeCount; nx = lNodeCount; endif

      xSize = lFaceGrid(2:end) - lFaceGrid(1:(end-1));
      assert(numel(xSize) <= ntx, sprintf('Number of nodes <%d> extracted from grid exceeds the number specified <%d>', numel(xSize), ntx));
      % Split the biggest node in two if there are fewer nodes than requested
      % To be optimized
      while numel(xSize) < ntx % The number of nodes is <numel(xSize)-1>
        % Double the size of xSize(1) and (end) if options "half node" is selected since their meaningful size is twice in that case
        [~,iw] = max(xSize + [xSize(1)*(xtype(1) == 'h'); ...
                              zeros(numel(xSize)-2,1); ...
                              xSize(end)*(xtype(2) == 'h')]);
        xSize(iw) /= 2;
        xSize = [xSize(1:iw); xSize(iw:end)];
      endwhile
    else % in the case of severals layers, it is a slightly more complex
      % We must checked the number of nodes in each layer
      xSize = lFaceGrid(2:end) - lFaceGrid(1:(end-1));
      xPos = cumsum(xSize) - 0.5*xSize;
      % Sort the nodes per layer
      xSizePerLayer = cell(Nlayer,1);

      for k=1:Nlayer
        xSizePerLayer{k} = xSize( (xPos >= Ltx_cum(k)/Ltx) && (xPos <= Ltx_cum(k+1)/Ltx) );
      endfor

      for k=1:Nlayer
        if isna(nx(k)), nx(k) = numel(xSizePerLayer{k}); endif

        xSizeLayer = xSizePerLayer{k};
        while numel(xSizeLayer) < nx(k)
          tmp = xSize;
          % Double the size of xSize(1) and (end) if options "half node" is selected since their meaningful size is twice in that case
          if (k==1) && (xtype(1) == 'h'), tmp(1) += xSize(1); endif;
          if (k==Nlayer) && (xtype(2) == 'h'), tmp(end) += xSize(end); endif;
          [~,iw] = max(xSize);
          xSize(iw) /= 2;
          xSize = [xSize(1:iw); xSize(iw:end)];
        endwhile

        xSizePerLayer{k} = xSizeLayer;
      endfor

      xSize = cell2mat(xSizePerLayer);
    endif

    xPos = cumsum(xSize) - 0.5*xSize;
    if (xtype(1) == 'h'), xPos(1) = 0; endif;
    if (xtype(2) == 'h'), xPos(end) = 1; endif;
  elseif size(xGrid, 2) > 1
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    assert(abs(sum(xGrid)-1) < 1E-12, 'Invalid mesh size vector. Sum should be 1');
    xSize = xGrid';
    xPos = cumsum(xSize) - 0.5*xSize;

    if (xtype(1) == 'h'), xPos(1) = 0; endif;
    if (xtype(2) == 'h'), xPos(end) = 1; endif;

    LayerPosCum = cumsum(Lx/Ltx);
    NodeSizeCum = cumsum(xSize);
    for i=1:Nlayer % Check that node sizes match with layer boundaries
      assert(min(abs(LayerPosCum(i) - NodeSizeCum)) < 1E-12, 'Node sizes do not match with layer sizes');
    endfor

    clear NodeSizeCum LayerPosCum;
  else % xGrid contains the node position
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    assert(abs(xGrid(end)-1) < 1E-12);
    assert(issorted(xGrid));
    xSize = [0; xGrid(:)];
    xSize = xSize(2:end) - xSize(1:(end-1));
    xPos = xGrid - 0.5*xSize;

    if (xtype(1) == 'h'), xPos(1) = 0; endif;
    if (xtype(2) == 'h'), xPos(end) = 1; endif;
  endif

  xDim = [cumsum([0; xSize(1:(end-1))]) , cumsum(xSize)];

endfunction
