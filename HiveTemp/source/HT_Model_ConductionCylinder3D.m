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
% Construit un modèle de conduction cylindrique 3D multilayer
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
% .radius = [1xNlayer] rayon des différentes épaisseur. Si le rayon intérieur est 0
%                      un cylindre plein est défini. C'est une condition spéciale
%                      puisqu'en l'état de l'implémentation, la conduction thermique
%                      n'est pas correctement modélisée dans cette configuration spéciale.
%                      Il sera donc supposé que tous les noeuds intérieurs sont
%                      à la même température. Ils sont ainsi tous réunis en un seul noeud.
%                      rGridType(1) doit être 'f'.
% .length = [1] longueur axial du cylindre
% .rhoC = [1xNlayer] capacité thermique de chaque couche
% .lambda = [3xNlayer] conductivités thermiques de chaque couche (r, theta, z)
% .Rc =  [1x(Nlayer-1)] résistance de contact
% .n = [NR1...xNR_Nlayer Ntheta Nz] discrétisation de chaque couche (nombre de noeuds)
% .rGrid = {[1 x NR1], [1 x NR2] ... [1 x NR_N]} spécifie la taille du maillage selon x  (normalized)
%          {[NR1 x 1], [NR2 x 1] ... [NR_N x 1]} spécifie la position des lignes du maillage  (normalized) (0 is implicit, not specified)
%          {Face, Face ... }
%          [Nlayer x 4] = [ radius_1, nNode_1, ratio_in1; ratio_out1; Specifies the radius, the number of node between this radius and the last one, and the condensation ratio
%                           radius_2, nNode_2, ratio_in2; ratio_out2;...]
% .rGridType = 'hh'/'hf'/'fh'/'ff' -> node type at XmXp position: 'h' for half node and 'f' for full node
% .zGrid = [1 x NY] spécifie la taille du maillage selon z   (normalized)
%          [NY x 1] spécifie la position des lignes du maillage  (normalized)
%          {cell/struct of faces} le maillage prendra en compte les noeuds des faces spécifiées.
% .zGridType = 'hh'/'hf'/'fh'/'ff' -> Only for extreme node (not between layers)
% .T0 = température initiale [1x1] ou [1xnNode]
% .axis = dim3x3 = [x y z] composantes des vecteurs x,y,z dans le repère global (défaut: Identity)
%                    Le cylindre est orienté vers Z+
% .globalPosition = [3x1] position of the wall (x,y,z)=(0,0,0) in
% .mergeTheta = [logical] Si <true>, on suppose qu'il n'y a pas de gradient de température suivant theta.
%                         Toutes les nodes sont fusionnées dans la matrice, toutefois, la géométrie est conservée.
% .base = [face] Face de référence utilisée pour extraire les données géométriques
%
% La structure <options> contient les champs suivant:
% .verbose
% .verboselevel
% .nodeNameModel : model name for nodes
% .skipFieldCheck : [true/false] prevent the parameter name check to be done
%
% .faces : retourne une structure contenant les champs suivant:
%     .top: liste des faces du dessus
%     .bottom: liste des faces du dessous
%     .inside: liste des faces internes
%     .outside: liste des faces externes
function [M, faces] = HT_Model_ConductionCylinder3D(name, params, options)
  HT_ImportConstants();

  lTicId = tic();
  lModelTypeName = 'model_conduction_cylinder_3D';

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
    lValidParamList = { 'radius', 'length', 'rhoC', 'lambda', 'Rc', 'n', 'material', ...
                       'nTheta', 'rGrid', 'rGridType', 'zGrid', 'zGridType', 'T0', 'axis', 'globalPosition', 'mergeTheta', 'base' };
    lParamList = fieldnames(params);
    lValidParam = cellfun(@(v) any(strcmp(v, lValidParamList)), lParamList);
    assert(all(lValidParam), sprintf('Invalid parameter name detected <%s>', strjoin(lParamList(~lValidParam), ',')));
    clear lParamList lValidParam lValidParamList;
  endif

  params = HT_CheckField(params, 'mergeTheta',       false, @(v) islogical(v));
  params = HT_CheckField(params, 'length',           [], {"exist", @(v) all(v > 0)});
  params = HT_CheckField(params, 'n',                NA(3,1));
  params = HT_CheckField(params, 'rGrid',            NA);
  params = HT_CheckField(params, 'axis',             NA(3,3), @(v) isempty(v) || all(size(v)==[3 3]));
  params = HT_CheckField(params, 'dim',              NA(2,1));
  params = HT_CheckField(params, 'globalPosition',   NA(3,1));
  params = HT_CheckField(params, 'radius',           []);
  params = HT_CheckField(params, 'base',             []);

  assert(isempty(params.base) || HT_CheckType(params.base, 'face'), 'Invalid parameter <base>');

  params = Int_ExtractParameterFromBase(params);

  if numel(params.radius) == 1, params.radius = [0 params.radius]; endif;
  if isrow(params.n), params.n = params.n(:); endif;

  lFullCylinder = params.radius(1) < 1E-10;
  Nlayer = numel(params.radius)-1;
  nrVec = params.n(1:(end-2));
  params = HT_CheckField(params, 'n', [], {"exist", @(v) (isfloat(v) && all((v > 0) & (round(v) == v)) && (numel(v) == Nlayer+2)) });

  % It allows a single face or a structure array of faces to be specified in xGrid field.
  % (not only cell array)
  % It avoid the use of double {{ }} that would be necessary
  % params.xGrid = {{ lMod_Underroof_faces{FaceZP} }}...
  if isfield(params, 'zGrid') && isstruct(params.zGrid)
    lFaceList = params.zGrid;
    params.zGrid = cell(numel(lFaceList), 1);
    for i=1:numel(lFaceList)
      params.zGrid{i} = lFaceList(i);
    endfor
  endif

  if isfield(params, 'rGrid') && isnumeric(params.rGrid)
    params.rGrid = { params.rGrid };
  endif

  params = HT_CheckField(params, 'radius',        [], {'exist', @(v) isfloat(v) && (all(v >= 0) || params.mergeTheta) && (numel(v) > 1)});
  params = HT_CheckField(params, 'length',        [], {'exist', @(v) isfloat(v) && (numel(v)==1) && (v > 0) });
  params = HT_CheckField(params, 'rGrid',         [], {@(v) isempty(v) || (iscell(v) && (numel(v) == Nlayer) && all(cellfun(@(x) isnumeric(x) || HT_CheckType(x, 'face') || isstruct(x),v)) ) });
  params = HT_CheckField(params, 'rGridType',     'hh', {@(v) any(strcmpi(v, {'hh', 'hf', 'fh', 'ff'}))});

  params = HT_CheckField(params, 'zGrid',         [], {@(v) isempty(v) || (numel(v) == params.n(end)) || iscell(v) || isstruct(v) });
  params = HT_CheckField(params, 'zGridType',     'hh', {@(v) any(strcmpi(v, {'hh', 'hf', 'fh', 'ff'}))});

  params = HT_CheckField(params, 'material',      [], @(v) all(HT_CheckType(params.material, 'material')) );
  params = Int_ExtractParameterFromMaterial(params);

  params = HT_CheckField(params, 'rhoC',        [], {'exist', @(v) isfloat(v) && all(v > 0) && any(numel(v) == [1 Nlayer])});
  params = HT_CheckField(params, 'lambda',      [], {'exist', @(v) isfloat(v) && all(v > 0) && any(rows(v) == [1 3]) && (columns(v) == Nlayer)});
  params = HT_CheckField(params, 'Rc',          zeros(Nlayer-1,1), @(v) isfloat(v) && all(v >= 0) && (any(numel(v) == [0 1 Nlayer-1])));
  params = HT_CheckField(params, 'T0',          NaN, @(v) isempty(v) || isscalar(v) || (numel(v) == sum(nrVec)));

  params = HT_CheckField(params, 'axis',        eye(3),   {@(v) all(size(v) == [3 3]) && all(vecnorm(v, Inf, 1) == 1.0)});
  assert(norm(cross(params.axis(:,1), params.axis(:,2)) - params.axis(:,3)) < 1E-15, 'Invalid <axis>');
  params = HT_CheckField(params, 'globalPosition', zeros(3,1),   {@(v) (numel(v) == 3)});

  assert((params.n(1) > 1)      || (params.rGridType(1) == 'f'), sprintf('rGridType <%s> is not compatible with 1 node only', params.rGridType));
  assert((params.n(end-2) > 1)  || (params.rGridType(2) == 'f'), sprintf('rGridType <%s> is not compatible with 1 node only', params.rGridType));
  assert((params.n(end) > 1)    || strcmpi(params.zGridType, 'ff'), sprintf('zGridType <%s> is not compatible with 1 node only', params.zGridType));

  % If internal radius is 0 (full cylinder)
##  assert(

  nr = sum(nrVec);
  ny = params.n(end-1);
  nz = params.n(end);
  nry = nr * ny;

  % Position, size, and boundaries of nodes
  [rPos, rSize, rDim]                 = Int_InitRadiusMeshSize(params.rGrid, params.n(1:(end-2)), params.radius, params.rGridType);
  [zPos, zSize, zDim]                 = Int_InitHeightMeshSize(params.zGrid, params.n(end), params.length, params.zGridType, params.axis(:,Z), params.globalPosition);

  lWallThickness = (params.radius(end)-params.radius(1));
  rDimTotal = [0; rDim(:,2) * lWallThickness] + params.radius(1);   % Position (not normalized) of each layer
  rDimTotalLog = log(rDimTotal(2:end)./rDimTotal(1:(end-1)));
  rDimTotalSqr = rDimTotal.^2;
  rNodePosTotal = params.radius(1) + rPos*lWallThickness;
  rNodeArea = 0.5 * sin(2*pi/ny) * (rDimTotalSqr(2:end) - rDimTotalSqr(1:(end-1)));
  %rNodeArea = sin(pi/ny) * (rDimTotal(1:(end-1)) + rDimTotal(2:end)) .* (rDimTotal(2:end) - rDimTotal(1:(end-1))) * cos(pi/ny); %rDimTotalSqr(2:end) - rDimTotalSqr(1:(end-1)));
  drVec = rSize * lWallThickness;
  dzVec = zSize * params.length;      % Absolute node size
  dzVecTotal = repelem(dzVec, nry);   % Absolute node size vector for all nodes
  if nz > 1
    dzNodeTotal = (zPos(2:end) - zPos(1:(end-1)) ) * params.length;
  endif

  M = HT_Model_Init(lModelTypeName, name, params);

  % Default value for nodeNameModel
##  if isempty(options.nodeNameModel), options.nodeNameModel = strcat(strrep(name, ' ', '_'), '.n%d'); endif;
  if isempty(options.nodeNameModel), options.nodeNameModel = strcat(strrep(name, ' ', '_'), '.n%d.%d.%d'); endif;

  % Count the number of occurence of %d
  idx = numel(strfind(options.nodeNameModel, '%d'));
  assert(any(idx == [1 3]), 'nodeNameModel must contains 1 or 3 occurences of <%d>');

  if idx == 1,
    NodeNameFunc = @(index, i,j,k) sprintf(options.nodeNameModel, index);
  elseif idx == 3,
    NodeNameFunc = @(index, i,j,k) sprintf(options.nodeNameModel, i,j,k);
  endif
  clear idx;

  if numel(params.rhoC) == 1, params.rhoC = repmat(params.rhoC, Nlayer, 1); endif;
  if numel(params.lambda) == 1, params.lambda = repmat(params.lambda, 1, Nlayer); endif;
  if rows(params.lambda) == 1, params.lambda = repmat(params.lambda, 3, 1); endif;
  assert(rows(params.lambda) == 3, 'Invalid matrix of conductivity');
  if numel(params.Rc) == 1, params.Rc = repmat(params.Rc, Nlayer-1, 1); endif;

  % G(node i <-> node j) = 2*pi/ny/lambda/log(rj/ri)
  % Calcul la demi resistance de chaque noeud
  RNodevec = [rDim(:,1) rPos rDim(:,2)] * lWallThickness + params.radius(1);
  RNodevec = [log(RNodevec(:,2)./RNodevec(:,1)) , log(RNodevec(:,3)./RNodevec(:,2))] ./ (2*pi*repelem(params.lambda(1,:)', nrVec, 1));

  % Take into account the subdivision of the disk
  % Ajoute les demi-resistances pour obtenir la résistance thermique noeud à noeud
  Rvec = (RNodevec(2:end,1) + RNodevec(1:(end-1),2))*ny;

  % Ajoute les résistances de contact si définies
  if ~isempty(params.Rc)
    lnrCum = cumsum(nrVec(1:(end-1)));
    Rvec(lnrCum) += params.Rc;

    clear lnrCum;
  endif

  % If mergeTheta is true, only one node is defined for the entire circle
  % The number of node is reduced by nTheta compared to the standard case
  if params.mergeTheta
    nNodes = nr * nz;

    % Vecteur des capacités et matrice des conductances
    % Calcul pour une part de gateau (voir cahier A p:17)
    % Multiplie par la capacité thermique de chaque couche
    Cvec = ny * rNodeArea .* repelem(params.rhoC(:), nrVec, 1);
    Cvec = repmat(Cvec, nz, 1);
    % Multiplie par la hauteur de chaque tranche horizontale
    Cvec .*= repelem(dzVec, nr, 1);

    % Construction des matrices
    % Conductance along r direction
    Gmat = spdiags( repmat([ny./Rvec; 0], nz, 1) .* repelem(dzVec, nr, 1), -1, nNodes, nNodes);

    % Conductance along z direction
    if nz > 1
      Gz = repelem(params.lambda(3,:)', nrVec, 1) .* (ny*rNodeArea);
      Gmat += spdiags(repmat(Gz, nz-1, 1) ./ repelem(dzNodeTotal, nr, 1),    -nr, nNodes, nNodes);
    endif

    Gmat += Gmat';          % Add the symetric part
    Gmat -= spdiags(sum(Gmat, 2), 0, nNodes, nNodes);   % Add the diagonal

    % Creation des noms des nodes
    M.nodes = cell(nNodes, 1);
    iNode = 0;
    for k=1:nz
      for i=1:nr
        iNode++;
        M.nodes{iNode} = NodeNameFunc(iNode, i,1, k);
      endfor
    endfor
  else
    nNodes = nry * nz;

    % Vecteur des capacités et matrice des conductances
    % Calcul pour une part de gateau (voir cahier A p:17)
    % Multiplie par la capacité thermique de chaque couche
    Cvec = rNodeArea .* repelem(params.rhoC(:), nrVec, 1);
    Cvec = repmat(Cvec, ny*nz, 1);
    % Multiplie par la hauteur de chaque tranche horizontale
    Cvec .*= repelem(dzVec, nry, 1);

    % Construction des matrices
    Gmat = sparse(nNodes,nNodes);

    % The following lines are much faster than the loop version
    Gmat += spdiags( repmat([1./Rvec; 0], ny*nz, 1) .* dzVecTotal, -1, nNodes, nNodes);

    % Condutance along theta direction
    Gy = repmat(repelem(params.lambda(2,:)', nrVec, 1) .* rDimTotalLog / (2*tan(pi/ny)), ny-1, 1);
    Gmat += spdiags( repmat([Gy; zeros(nr,1)], nz, 1) .* dzVecTotal , -nr, nNodes, nNodes);
    Gmat += spdiags( repmat([Gy(1:nr); zeros(nry-nr,1)], nz, 1) .* dzVecTotal , -(nry-nr), nNodes, nNodes);

    % Conductance along z direction
    if nz > 1
      Gz = repelem(params.lambda(3,:)', nrVec, 1) .* rNodeArea;
      Gmat += spdiags(repmat(Gz, ny*(nz-1), 1) ./ repelem(dzNodeTotal, nry, 1),    -nry, nNodes, nNodes);
    endif

    Gmat += Gmat';          % Add the symetric part
    Gmat -= spdiags(sum(Gmat, 2), 0, nNodes, nNodes);   % Add the diagonal

    % Creation des noms des nodes
    M.nodes = cell(nNodes, 1);
    iNode = 0;
    for k=1:nz
      for j=1:ny
        for i=1:nr
          iNode++;
          M.nodes{iNode} = NodeNameFunc(iNode, i,j, k);
        endfor
      endfor
    endfor

    % Full cylinder means inner radius is 0. In that case, inner nodes are
    % merged together to allow heat to cross the center of the cylinder
    % Otherwise, the thermal resistance become infinite.
    if lFullCylinder
      % Replace internal nodes names but keeps the cell array untouched
      % to simplify the face construction (see below). The node cell array
      % is changed at the end of this function
      for k=1:nz
        lIndexRef = (k-1)*nr*ny+1;
        lCenterNodes = lIndexRef+(nr:nr:(nry-nr));
        M.nodes(lCenterNodes) = M.nodes{lIndexRef};
      endfor

      % Matrices are changed at the end of the function
    endif

##  Rvec = [0; t; 0];
##  Gr = spdiags(1./Rvec(1:(end-1)),                                                      1, nr, nr);
##  Gy = spdiags(repelem(params.lambda(2,:)(:), nrVec) .* rPosTotalLog / (2*tan(pi/ny)),  0, nr, nr);
##  Gz = spdiags(repelem(params.lambda(3,:)(:), nrVec) .* rNodeArea,                      0, nr, nr);
##
##  GmatChk = sparse(nNodes,nNodes);
##  % Fill only the top right part of the matrix
##  for k=1:nz
##    planeRef = (k-1)*nry;
##
##    for j=1:ny
##      nodeRef = (k-1)*nry + (j-1)*nr;
##      nodeSet = nodeRef + (1:nr);
##      GmatChk(nodeSet, nodeSet) = Gr * dzVec(k);
##
##      % Connect this part to the next
##      if j < ny
##        GmatChk(nodeSet, nr+nodeSet) = Gy * dzVec(k);
##      else
##        GmatChk(nodeSet, planeRef+(1:nr)) = Gy * dzVec(k);
##      endif
##
##      if k < nz
##        GmatChk(nodeSet, nry+nodeSet) = Gz / dzNodeTotal(k);
##      endif
##    endfor
##  endfor
##
##  GmatChk += GmatChk';          % Add the symetric part
##  GmatChk -= spdiags(sum(GmatChk, 2), 0, nNodes, nNodes);   % Add the diagonal
##  X = Gmat-GmatChk;
##  disp(max(max(abs(full(X)))));
  endif

  M.G = Gmat;
  M.C = Cvec;
  M.T0 = params.T0;

  if isempty(M.T0), M.T0 = NaN(nNodes,1); endif
  if isscalar(M.T0), M.T0 = repmat(M.T0, nNodes, 1); endif;
  assert(numel(M.T0) == nNodes, 'Invalid dimension of <.T0>');

  M.axis = params.axis;
  M.globalPosition = params.globalPosition(:);
  clear nodes Gmat iNode Cvec Gy Rvec;

  if nargout > 1
    faces = struct('top', [], 'bottom', [], 'inside', [], 'outside', []);
    lRadiusCorrection = cos(pi/ny);
    lRadius = params.radius(end);
    lRealRadius = params.radius(end)*lRadiusCorrection;

    nyNodes = ny;
    if params.mergeTheta, nyNodes = 1; endif % nyNodes is used to select node names in M.nodes.

    % Top faces
    faces.top = HT_Face_Init( strcat(name, '.top'), ...
                              'size', repmat(2*lRealRadius, 2, 1), ...
                              'nodes', {}, ...
                              'globalPosition', M.globalPosition + params.length * M.axis(:,3) - lRealRadius*sum(M.axis(:,1:2), 2), ...
                              'axis', [M.axis(:,1) M.axis(:,2)], ... xy;
                              'norm', M.axis(:,3), ...  Top face oriented toward Z+
                              'material', params.material, ...
                              'model', name,
                              'metadata', HT_Face_InitMetaData( 'type', 'circle', ...
                                                                'uSize', params.radius, ...
                                                                'vSize', pi/ny + [0 2*pi], ...
                                                                'uGrid', [0; rDim(:,2)], ...
                                                                'vGrid', int32(ny), ...
                                                                'uNode', nrVec, ...
                                                                'uNodePos', rPos));

    lVertexUnit = [cos((2*pi/ny)*((1:ny)' - 0.5)) sin((2*pi/ny)*((1:ny)' - 0.5))];
##    lVertexUnitLoop = [cos((2*pi/ny)*((0:ny)' - 0.5)) sin((2*pi/ny)*((0:ny)' - 0.5))];
    lVertexUnitLoop = [lVertexUnit(end,:) ; lVertexUnit];

    lNodeUnit = [cos((2*pi/ny)*((1:ny)' - 1)) sin((2*pi/ny)*((1:ny)' - 1))];
##    lNodeUnitLoop = [cos((2*pi/ny)*((0:ny)' - 1)) sin((2*pi/ny)*((0:ny)' - 1))];
    lNodeUnitLoop = [lNodeUnit(end,:) ; lNodeUnit];

    faces.top.pos = repmat(rNodePosTotal * lRadiusCorrection, ny, 1) .* repelem(lNodeUnit, nr ,1);
    faces.top.pos += lRealRadius;
    faces.top.pos ./= repmat(2*lRealRadius, 1, 2);
    faces.top.vertex = repmat(rDimTotal, ny, 1) .* repelem(lVertexUnit, (nr+1), 1);
    faces.top.vertex += lRealRadius;
    faces.top.vertex ./= repmat(2*lRealRadius, 1, 2);
    faces.top.dims = int32( (1:(nr*ny))' + repelem((0:(ny-1))',nr,1) + ...
                                [ 0, -nr-1, -nr, 1]); %Vertices are defined the direct order to match the normal vector
    faces.top.dims(1:nr,2:3) += (nr+1)*ny;
    faces.top.edges = int32([(1:ny)*(nr+1)    , nr+1 ; ... External lines
                             (1:ny)*(nr+1)-nr , 1 ]); ... Internal lines
    faces.top.polygons = int32( 1 + (0:(ny-1))'*(nr+1) + ...
                                [ 0, -nr-1, -1, nr]); ... %Vertices are defined the direct order to match the normal vector
    faces.top.polygons(1,2:3) += (nr+1)*ny;
    faces.top.r = zeros(nry, 1);

    % Compute the half thermal resistance of the last layer of nodes
    if params.zGridType(2) == 'f'
      faces.top.r = 0.5*dzVec(end) ./ repmat(repelem(params.lambda(3,:)', nrVec, 1), ny, 1); %./ repmat(Gz, ny, 1);
    endif


    faces.bottom = HT_Face_Init( strcat(name, '.bottom'), ...
                              'size', repmat(2*lRealRadius, 2, 1), ...
                              'nodes', {}, ...
                              'globalPosition', M.globalPosition - lRealRadius*(M.axis(:,1) - M.axis(:,2)), ...
                              'axis', [M.axis(:,1) -M.axis(:,2)], ... (x,-y);
                              'norm', -M.axis(:,3), ...  Top face oriented toward Z+
                              'material', params.material, ...
                              'model', name, ...
                              'metadata', HT_Face_InitMetaData( 'type', 'circle', ...
                                                                'uSize', params.radius, ...
                                                                'vSize', pi/ny + [0 2*pi], ...
                                                                'uGrid', [0; rDim(:,2)], ...
                                                                'vGrid', int32(ny), ...
                                                                'uNode', nrVec, ...
                                                                'uNodePos', rPos));

    faces.bottom.pos = faces.top.pos;
    faces.bottom.pos(:,2) = 1 - faces.bottom.pos(:,2);
    faces.bottom.vertex = faces.top.vertex;
    faces.bottom.vertex(:,2) = 1 - faces.bottom.vertex(:,2);
    faces.bottom.dims = faces.top.dims(:,[1 4 3 2]); % Invert the orientation of the quad
    faces.bottom.edges = faces.top.edges;
    faces.bottom.polygons = faces.top.polygons(:,[1 4 3 2]);
    faces.bottom.r = zeros(nry, 1);

    if params.zGridType(1) == 'f'
      faces.bottom.r = 0.5*dzVec(1) ./ repmat(repelem(params.lambda(3,:)', nrVec, 1), ny, 1); % ./ repmat(Gz, ny, 1);
    endif

    if params.mergeTheta
      faces.top.nodes = repmat(M.nodes((nr*(nz-1)+1):nNodes), ny, 1); % Duplicate node names since different
                                                                      % quads refer to the same node

      faces.bottom.nodes = repmat(M.nodes(1:nr), ny, 1); % Duplicate node names since different
                                                         % quads refer to the same node
    else
      faces.top.nodes = M.nodes((nry*(nz-1)+1):nNodes);
      faces.bottom.nodes = M.nodes(1:nry);
    endif

    if Nlayer > 1
      faces.top.materialIndex = repmat( repelem((1:Nlayer)', nrVec) , ny, 1);
      faces.bottom.materialIndex = repmat( repelem((1:Nlayer)', nrVec) , ny, 1);
    endif

    % Inside faces
    if params.radius(1) > 1E-10
      faces.inside = HT_Face_Init( cellstr(strcat(name, '.in', strrep(int2str((1:ny)'), ' ', ''))), ...
                                'size', [2*params.radius(1)*sin(pi/ny) params.length], ...
                                'n', [1 nz], ...
                                'nodes', M.nodes(1:nr:end), ...
                                'material', params.material(1), ...
                                'r', 0, ...
                                'model', name);
    endif

    faces.outside = HT_Face_Init( cellstr(strcat(name, '.out', strrep(int2str((1:ny)'), ' ', ''))), ...
                              'size', [2*params.radius(end)*sin(pi/ny) params.length], ...
                              'n', [1 nz], ...
                              'nodes', M.nodes(1:nr:end), ...
                              'material', params.material(end), ...
                              'r', 0, ...
                              'model', name);

##      figure(1); clf; hold on;

    for i=1:ny

      faces.outside(i).nodes = M.nodes(nr + (i-1)*nr*(~params.mergeTheta) +(0:(nz-1))*(nyNodes*nr));
      faces.outside(i).globalPosition = M.globalPosition + M.axis*[lVertexUnitLoop(i,:)' * params.radius(end) ; 0];
      t = (lVertexUnitLoop(i+1,:) - lVertexUnitLoop(i,:))'; % Local u axis
      t /= vecnorm(t);
      faces.outside(i).axis = [M.axis * [t; 0] , M.axis(:,Z) ]; ;
      faces.outside(i).norm = M.axis * [lNodeUnit(i,:)'; 0];
      faces.outside(i).pos = [repmat(0.5, nz, 1), zPos];
      faces.outside(i).dims = [zeros(nz,1) ones(nz, 1) zDim];

      % Add half node radius
      if params.rGridType(2) == 'f'
        % Add the half node thermal resistance
        % Convert resistance to normalized resistance (K.m²/W)
        % Size along z axis was to introduced in RNodevec. Therefor, it does not appear here.
        faces.outside(i).r = repmat(RNodevec(end,2), nz, 1)*faces.outside(i).size(U);
      endif

##      C = HT_Face_GetAbsoluteData(faces.outside(i), 'corners');
##
##      h = line(   "xdata", [C(1,:), C(1,1)], ...
##                  "ydata", [C(2,:), C(2,1)], ...
##                  "zdata", [C(3,:), C(3,1)], ...
##                  "linewidth", 1, ...
##                  "color", [0 0 0]);

      if params.radius(1) > 1E-10
        faces.inside(i).nodes = M.nodes(1+(i-1)*nr*(~params.mergeTheta) +(0:(nz-1))*(nyNodes*nr));
        faces.inside(i).globalPosition = M.globalPosition + M.axis*[lVertexUnit(i,:)' * params.radius(1); 0];
        faces.inside(i).axis = [-faces.outside(i).axis(:,1), faces.outside(i).axis(:,2)]; % Global Z vector is local v vector
        faces.inside(i).norm = -faces.outside(i).norm;
        faces.inside(i).pos = faces.outside(i).pos;
        faces.inside(i).dims = faces.outside(i).dims;

        % Add half node radius
        if params.rGridType(1) == 'f'
          % Add the half node thermal resistance
          % Convert resistance to normalized resistance (K.m²/W)
          % Size along z axis was to introduced in RNodevec. Therefor, it does not appear here.
          faces.inside(i).r = repmat(RNodevec(1,1), nz, 1)*faces.inside(i).size(U);
        endif
      endif

##      C = HT_Face_GetAbsoluteData(faces.inside(i), 'corners');
##
##      h = line(   "xdata", [C(1,:), C(1,1)], ...
##                  "ydata", [C(2,:), C(2,1)], ...
##                  "zdata", [C(3,:), C(3,1)], ...
##                  "linewidth", 1, ...
##                  "color", [0 0 0]);
    endfor

##  figure(1);
##  clf;
##  X=[faces.inside(1).globalPosition]; plot3(X(1,:)', X(2,:)', X(3,:)', '+r', 'markersize', 10);
##  X=[faces.outside(1).globalPosition]; plot3(X(1,:)', X(2,:)', X(3,:)', '+r', 'markersize', 10);
##  xlabel('X');
##  ylabel('Y');

endif

  % Option full cylinder ?
  % Faces are not built. Duplicate nodes in the cell array M.nodes are now removed
  % and matrices are changed accordingly.
  if ~params.mergeTheta && lFullCylinder

    lMergeNodes = cell(nz, 1);
    lNewNames = cell(nz, 1);
    for k=1:nz
      lMergeNodes{k} = (k-1)*nr*ny+1 + (0:nr:(nry-1));
      lNewNames{k} = M.nodes{lMergeNodes{k}(1)};
    endfor

    M = HT_Model_MergeNodes(M, lMergeNodes, lNewNames);

##    lToBeRemoved = false(nNodes, 1);
##
##    for k=1:nz
##      lRefIndex = (k-1)*nr*ny+1;
##      lCenterNodes = false(nNodes, 1);
##      lCenterNodes(lRefIndex+(nr:nr:(nry-1))) = true;
##
##      lMergeLine =    sum(M.G(lCenterNodes,:), 1);
##      lMergeColumn =    sum(M.G(:,lCenterNodes), 2);
##      M.G(lRefIndex,:) += lMergeLine;
##      M.G(:,lRefIndex) += lMergeLine';
##      M.G(lRefIndex,lRefIndex) = 0;
##      M.G(lRefIndex,lRefIndex) = -sum(M.G(lRefIndex,~lCenterNodes));
##
##      lToBeRemoved = lToBeRemoved | lCenterNodes;
##      M.G(lCenterNodes,:) = 0;
##      M.G(:, lCenterNodes) = 0;
##    endfor
##
##    M.G(lToBeRemoved,:) = [];
##    M.G(:,lToBeRemoved) = [];
  endif

  M.timer.buildingTime += toc(lTicId);

  if options.verbose, disp(sprintf('Building models <%s> done in %.2f ms', name, 1000*M.timer.buildingTime)); endif; % Start timer (tic)

endfunction


##% nx [Nlayer x 1] number of nodes per layer
##% Lx [Nlayer x 1] length of each layer
##% xtype = 'xx' with x=h/f for node type at boundaries
##% vDir [dim 3x1] direction in global axis
##% vPos [dim 3x1] global position
##% Returns:
##% xPos = node position (normalized)
##% xSize = node size (normalized)
##% xDim = [Nx1] [xstart xend; ...] with N=sum(nx)
##% nx = [Nx1] the number of nodes for each layer
function [rPos rSize rDim] = Int_InitRadiusMeshSize(rGrid, nr, Lr, rType)
  global HT_VAR_EPSILON_U;

  assert(numel(nr)+1 == numel(Lr), 'Invalid size of node count vector. Must be equal to the number of layers.');
  assert(numel(rGrid) == numel(nr), 'Invalid size of node count vector. Must be equal to the number of layers.');

  Nlayer = numel(nr);
  Ltr = Lr(end) - Lr(1);
  lLayerThickVec = Lr(2:end) - Lr(1:(end-1));
  ntr = sum(nr);
  ntrVec = [0; cumsum(nr)];

  fdr = @(m) 1 / (nr(m) - 0.5*((m==1) && (rType(1)=='h')) - 0.5*((m==Nlayer) && (rType(2)=='h')));
##  fdr = @(m) lLayerThickVec(m) / Ltr / (nr(m) - 0.5*((m==1) && (rType(1)=='h')) - 0.5*((m==Nlayer) && (rType(2)=='h')));

  rPos = NA(ntr+1, 1);
  rSize = NA(ntr, 1);

  % Used with condense parameter (to reduce node size around boundaries)
  f1 = @(x,k) x.^k;
  f2 = @(x,k) 1-(1-x).^k;

  for i=1:Nlayer
    lrGrid = rGrid{i};
    lrSet = (ntrVec(i)+1):ntrVec(i+1);

    if isempty(lrGrid) % Not defined by user ? Use default regular grid
      lrGrid = repmat(fdr(i), 1, nr(i));
      if i==1 && (rType(1)=='h'), lrGrid(1) *= 0.5; endif; % Half node ?
      if i == Nlayer && (rType(2)=='h'), lrGrid(end) *= 0.5; endif; % Half node ?
    elseif isstruct(lrGrid) && ~HT_CheckType(lrGrid, 'face')
      lrGridInfo = lrGrid;

      t = repmat(fdr(i), 1, nr(i));
      if i==1 && (rType(1)=='h'), t(1) *= 0.5; endif; % Half node ?
      if i == Nlayer && (rType(2)=='h'), t(end) *= 0.5; endif; % Half node ?
      t = [0, cumsum(t)];
      t = t(:);

      lrGridInfo = HT_CheckField(lrGridInfo, 'e1',             1, {@(v) isscalar(v) && v >= 0});
      lrGridInfo = HT_CheckField(lrGridInfo, 'e2',             1, {@(v) isscalar(v) && v >= 0});

      if lrGridInfo.e2 == 0
        t = f1(t, lrGridInfo.e1);
      elseif lrGridInfo.e1 == 0
        t = f2(t, lrGridInfo.e2);
      else
        m=2;
        t = (f1(t,lrGridInfo.e1).*(1-t).^m + f2(t,lrGridInfo.e2).*t.^m)./((1-t).^m+t.^m);
      endif

      lrGrid = round(t/HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;

##      lrGrid = t(2:end)-t(1:(end-1));
    endif

    if HT_CheckType(lrGrid, 'face')
      M = lrGrid.metadata;
      % Extract information from the face metadata
      assert(strcmpi(M.type, 'circle'), 'Invalid face type. Must be a circle');
      FrDimTotal = M.uSize(1) + (M.uSize(end)-M.uSize(1)) * M.uGrid;
##      FrDimTotal = round(FrDimTotal /HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS;

      % Select the grid corresponding to this layer
      FrDimTotal((FrDimTotal < Lr(i)) | (FrDimTotal > Lr(i+1))) = [];
      FrDimTotal = ([Lr(i); FrDimTotal; Lr(i+1)] - Lr(1))/Ltr;
      FrDimTotal = unique(round(FrDimTotal /HT_VAR_EPSILON_U)*HT_VAR_EPSILON_U);

      assert(nr(i) == numel(FrDimTotal)-1, sprintf('Invalid node count for layer %d. Specified <%d> but <%d> are defined in face <%s>', i, nr(i), numel(FrDimTotal)-1, lrGrid.name));
      rSize(lrSet) = FrDimTotal(2:end) - FrDimTotal(1:(end-1));
      % TODO: node position is not properly extracted and may not match the position
      % in the base face, face the error is neglected for the moment.

    elseif rows(lrGrid) == 1 % Mesh size defined
      if abs(sum(lrGrid) -1) < 1E-10
        lrGrid *= lLayerThickVec(i) / Ltr;  % If the normalized node sizes were normalized for each layer, it is  rescaled.
      else
        assert(sum(lrGrid) < 1.0, sprintf('Invalid node size of layer %d', i));
      endif

      rSize(lrSet) = lrGrid;
##      rPos(lrSet(1)) = (Lr(i)-Lr(1))/Ltr;
##      rPos(lrSet+1) = rPos(lrSet(1)) + cumsum(lrGrid(:));
    elseif columns(lrGrid) == 1 % Node boundaries defined
      assert(numel(lrGrid) == nr(i)+1, sprintf('Invalid size of mesh position vector. Should be %d, but %d is specified', nr(i)+1, numel(lrGrid)));

      if abs(lrGrid(end)-1.0) < 1E-10 % If the normalized position intervalle is [0;1] for each layer, rescale it...
        assert((i==1) || (abs(lrGrid(1)) < 1E-10), sprintf('Invalid mesh position vector. First element should be 0 at layer %d', i));

        lrGrid *= lLayerThickVec(i) / Ltr;
        lrGrid += (Lr(i)-Lr(1))/Ltr;
      else
        assert(i < Nlayer, sprintf('Invalid mesh position vector. Last element should be 1 at layer %d but %.2e specified', i, lrGrid(end)));
        assert(isempty(rGrid{i+1}) || (abs(lrGrid(end) - rGrid{i+1}(1)) < 1E-10), sprintf('Invalid mesh position between layer %d and %d', i, i+1));
      endif

      rSize(lrSet) = lrGrid(2:end) - lrGrid(1:(end-1));
##      rPos(lrSet(1)) = lrGrid(1);
##      rPos(lrSet+1) = lrGrid(2:end);
    endif

  endfor

  rSizeCum =  cumsum(rSize);
  rPos = rSizeCum - 0.5*rSize;

  if (rType(1)=='h'), rPos(1) = 0; endif; % Half node ?
  if (rType(2)=='h'), rPos(end) = 1; endif; % Half node ?

  assert(abs(rSizeCum(end)-1) < 1E-10, 'Internal error');

  rDim = [[0; rSizeCum(1:(end-1))] , rSizeCum];

endfunction

function [xPos xSize xDim, nx] = Int_InitHeightMeshSize(xGrid, nx, Lx, xType, vDir, vPos)
  global HT_VAR_EPSILON_POS;
  global HT_VAR_EPSILON_U;

  assert(numel(nx) == numel(Lx), 'Invalid size of node count vector. Must be equal to the number of layers.');
  Nlayer = numel(nx);
  Ltx = sum(Lx);
  Ltx_cum = cumsum([0; Lx]);
  ntx = sum(nx);

  xPos = [];
  xSize = [];

  if isempty(xGrid)
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    fdx = @(m) Lx(m) / Ltx / (nx(m) - 0.5*((m==1) && (xType(1)=='h')) - 0.5*((m==Nlayer) && (xType(2)=='h')));

    for i=1:Nlayer
      ldx = fdx(i);
      lxSize = repmat(ldx, nx(i), 1);
      if i == 1 && (xType(1)=='h'), lxSize(1) *= 0.5; endif; % Half node ?
      if i == Nlayer && (xType(2)=='h'), lxSize(end) *= 0.5; endif; % Half node ?

      lxPos = sum(Lx(1:(i-1))/Ltx) + ((1:nx(i))'-0.5)*ldx - (0.5*ldx * ((xType(1)=='h') && (i==1)));

      xSize = [xSize ; lxSize];
      xPos = [xPos ; lxPos];
    endfor

    assert(abs(xPos(end)-1+0.5*ldx*(xType(1)!='h')) < 1E-12);
    assert(abs(sum(xSize)-1) < 1E-12);
    clear ldx lxSize lxPos;
##  elseif iscell(xGrid) && all(cellfun(@(v) isnumeric(v), xGrid))
##    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
##    % Converts coordinates in xGrid vector
##    if size(xGrid, 2) > 1, xGrid = xGrid'; endif;
##
##    xGridpos = cell2mat(xGrid);
##    xGridn = xGridpos(2:2:end);
##    xGridpos = xGridpos(1:2:end);
##    assert(all(xGridn > 0), 'Invalid node count in xyzgrid vector');
##    assert(all(xGridpos > 0), 'Invalid node position in xyzgrid vector');
##
##    lOneInd = find(abs(xGridpos-1) < 1E-14);
##    if numel(lOneInd) == Nlayer % Ex: { 0.0, 4, 0.5, 4, 1.0, 10, 1.0, 8 } % 1.0 indicates next layer
##      for i=1:numel(lOneInd)
##        if i == 1
##          xGridpos(1:lOneInd(i)) *= Lx(i) / Ltx;
##          assert(sum(xGridn(1:lOneInd(i))) == nx(i), 'Wrong node count');
##        else
##          xGridpos((lOneInd(i-1)+1):lOneInd(i)) = xGridpos(lOneInd(i-1)) + ...
##                                                  xGridpos((lOneInd(i-1)+1):lOneInd(i)) * Lx(i)/Ltx;
##          assert(sum(xGridn((lOneInd(i-1)+1):lOneInd(i))) == nx(i), 'Wrong node count');
##        endif
##
##      endfor
##
##      clear lOneInd;
##    else
##      assert(numel(lOneInd) == 1, 'Pos 1.0 should appear once (or nlayer)');
##      xGridpos = xGrid(1:2:end);
##      xGridn = xGrid(2:2:end);
##    endif
##
##    assert(issorted(xGridpos));
##
##    % xGrid contient les positions relatives et les nombres de noeuds
##    fdx = @(L, m, mtotal) L / (xGridn(m) - 0.5*((m==1) && (xType(1)=='h')) - 0.5*((m==mtotal) && (xType(2)=='h')));
##
##    % Cas particulier pour le premier intervalle
##    nPart = numel(xGridpos);
##    ldx = fdx(xGridpos(1), 1, nPart);
##    lxSize = repmat(ldx, xGridn(1), 1);
##    if xType(1)=='h', lxSize(1) *= 0.5; endif; % Half node ?
##    if ((nPart==1) && (xType(2)=='h')), lxSize(end) *= 0.5; endif; % Half node ?
##    xSize = lxSize;
##    xPos = ((1:xGridn(1))'-0.5)*ldx - (lxSize(1) * (xType(1)=='h'));
##
##    for i=2:nPart
##      ldx = fdx(xGridpos(i)-xGridpos(i-1), i, nPart);
##      lxSize = repmat(ldx, xGridn(i), 1);
##      if i == nPart && (xType(2)=='h'), lxSize(end) *= 0.5; endif; % Half node ?
##
##      lxPos = xGridpos(i-1) + ((1:xGridn(i))'-0.5)*ldx;
##
##      xSize = [xSize ; lxSize];
##      xPos = [xPos ; lxPos];
##    endfor
##
##    assert(abs(xPos(end)-1+0.5*ldx*(xType(1)!='h')) < 1E-12);
##    assert(abs(sum(xSize) -1) < 1E-12);
##    clear ldx lxSize lxPos nPart;
  elseif iscell(xGrid) && all(cellfun(@(v) HT_CheckType(v, 'face'), xGrid))
% nx [Nlayer x 1] number of nodes per layer
% Lx [Nlayer x 1] length of each layer
% xType = 'xx' with x=h/f for node type at boundaries
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
      lGridPos = HT_Face_GetAbsoluteData(F, 'xSet');
      % And convert the row vector to column vector
      lFaceGrid = [lFaceGrid; lGridPos];
    endfor

    lFaceGrid /= Ltx; % normalized

    % Remove invalid values (if face is higher than current model)
    lFaceGrid((lFaceGrid < 0) | (lFaceGrid > 1)) = [];
    % Now add the layers
    lFaceGrid = [lFaceGrid; 0; Lx/Ltx];
    lFaceGrid = floor(lFaceGrid / HT_VAR_EPSILON_U + 0.5)*HT_VAR_EPSILON_U;
    lFaceGrid = unique(lFaceGrid); % Unique and sorted

    % The final number of nodes is either the number specified by the user
    % or the number infered by the list of faces.
    lNodeCount = max(ntx, numel(lFaceGrid)-1);

    % If the number of nodes was not specified
    % <xN> is set and returned to the main function
    if isna(ntx), ntx = lNodeCount; nx = lNodeCount; endif

    xSize = lFaceGrid(2:end) - lFaceGrid(1:(end-1));
    assert(numel(xSize) <= ntx, sprintf('Number of nodes <%d> extracted from faces exceeds the number specified in params structure <%d>', numel(xSize), ntx));
    % Split the biggest node in two if there are fewer nodes than requested
    % To be optimized
    while numel(xSize) < ntx % The number of nodes is <numel(xSize)-1>
      % Double the size of xSize(1) and (end) if options "half node" is selected since their meaningful size is twice in that case
      [~,iw] = max(xSize + [xSize(1)*(xType(1) == 'h'); ...
                            zeros(numel(xSize)-2,1); ...
                            xSize(end)*(xType(2) == 'h')]);
      xSize(iw) /= 2;
      xSize = [xSize(1:iw); xSize(iw:end)];
    endwhile

    xPos = cumsum(xSize) - 0.5*xSize;
    if (xType(1) == 'h'), xPos(1) = 0; endif;
    if (xType(2) == 'h'), xPos(end) = 1; endif;
  elseif size(xGrid, 2) > 1
    assert(~any(isna(nx)), 'The number of nodes is not fully defined');
    assert(abs(sum(xGrid)-1) < 1E-12, 'Invalid mesh size vector. Sum should be 1');
    xSize = xGrid';
    xPos = cumsum(xSize) - 0.5*xSize;

    if (xType(1) == 'h'), xPos(1) = 0; endif;
    if (xType(2) == 'h'), xPos(end) = 1; endif;

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

    if (xType(1) == 'h'), xPos(1) = 0; endif;
    if (xType(2) == 'h'), xPos(end) = 1; endif;
  endif

  xDim = [[0; cumsum(xSize(1:(end-1)))] , cumsum(xSize)];

endfunction

function params = Int_ExtractParameterFromBase(params)
  F = params.base;

  if any(isna(params.globalPosition))
    lFaceCenter = HT_Face_GetAbsoluteData(F, 'center');
    s = isna(params.globalPosition);
    params.globalPosition(s) = lFaceCenter(s);
  endif

  if any(isna(params.axis))
    lFaceAxis = [F.axis, F.norm];
    s = isna(params.axis);
    params.axis(s) = lFaceAxis(s);
  endif
endfunction

function params = Int_ExtractParameterFromMaterial(params)
  if isempty(params.material), return; endif
  lNlayer = numel(params.radius)-1;
  assert(numel(params.material) == lNlayer, 'Invalid parameter <material>. Size does not match the layer count.');

  if iscell(params.material)
    params.rhoC = cellfun(@(v) HT_Material_GetRhoC(v), params.material);
    if iscolumn(params.rhoC) params.rhoC = params.rhoC'; endif

    params.lambda = cellfun(@(v) HT_Material_GetLambda(v, []), params.material);
  elseif isstruct(params.material)
    params.rhoC = arrayfun(@(v) HT_Material_GetRhoC(v), params.material);
    params.lambda = cell2mat(arrayfun(@(v) HT_Material_GetLambda(v), params.material(:)', 'Uniformoutput', false));
  else
    error('Invalid material');
  endif
endfunction

