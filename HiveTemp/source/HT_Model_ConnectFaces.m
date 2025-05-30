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
% This function defines a thermal connection between two faces <F1> and <F2>
% It returns a thermal model M with name <name>
%
% Input parameters:
% .name: [string] the new model name
% .F1: [face object] the first face to connect
% .F2: [face object] the second face to connect
% .params: [parameters structure] contains all necessary parameters
%     -> .g: [scalar/column vector/function handle] the conductance (W/K/m²)
%            Specifying a column vector allows to define a space varying conductance. In that
%            case, the number of elements in g must be equal to the number of nodes and
%            must match the order of F1 object.
%            Specifying a function handle allows to define a time varying conductance. Moreover
%            it can be space varying as well. It can be defined as @(t,i) or @(T,t,i)
%            "t" is the current time, and "i" the current time index
%            If input argument "T" is used, it means the conductance is non linear
%            "T" will be provided during computations. For a Node to Node connection,
%            "T" will be dimension 2x1 (the temperature for each node will be provided)
%            For nNodes nodes, "T" will be of dimension "nNodes+1" x 1
%     -> .T0: [scalar] used to specify initial temperature for all nodes related to this
%            connection (F1.nodes and F2.nodes)
%     -> .uv1: [vector 1x4] [ustart, uend, vstart, vend] is used to specify a fraction
%            of face 1.
%     -> .uv2: [vector 1x4] [ustart, uend, vstart, vend] is used to specify a fraction
%            of face 2.
%     -> .uv1axis: [vector 3x2] axis in general axix system of u et v vector associated
%            to parameter <.uv1>
%     -> .uv2axis: [vector 3x2] axis in general axix system of u et v vector associated
%            to parameter <.uv2>
%     -> .intersect: [logical] [false] If true, the intersection of both faces is found first
%     -> .epsilon: [scalar] default is 1E-10. Value used to compare vertices
% .options: [parameters structure] Structure used to specify options.
%     -> .verbose: [boolean] default is false. If true, display information.
%     -> .checkGeometry: [boolean] default is true. If true, geometry checks is performed
%            during connection to check that face indeed connected and have the same orientation.
%     -> .strictMatch: [boolean] default is true. If true, an error is generated
%            when nodes of two different faces does not perfectly match, i.e. they
%            don't have exactly the same size. Unmatch node connexions could
%            introduce some bias.
%
% Output parameters:
% -> M: [model object] Returns the model associated to this connection
%

##-*- texinfo -*-
##@deftypefn {} {@var{M} =} HT_Model_ConnectFaces (@var{name}, @var{F1}, @var{F2}, @var{options})
##@cindex model face
##This function connects two faces and returns the corresponding thermal model.
##
##like @code{sample of code} and variables should be marked
##as @var{variable}.
##@seealso{HT_Model_Init, HT_Model_Connect}
##@end deftypefn
##
## Author: Emmanuel Ruffio <emmanuel.ruffio@alt-rd.com>
##
function M = HT_Model_ConnectFaces(name, F1, F2, params, options)
  HT_ImportConstants();

  lTicId = tic();
  lModelTypeName = 'model_face_connection';

  lParams = params;
  lParams.face1 = F1;
  lParams.face2 = F2;
  M = HT_Model_Init(lModelTypeName, name, lParams);
  clear lParams;

  assert(nargin >= 3, 'There must be 4 input arguments');

  if nargin < 5, options = struct(); endif;
  if nargin < 4, params = struct(); endif;

  if isempty(params), params = struct(); endif;
  if isempty(options), options = struct(); endif;

  options = HT_CheckField(options, 'verbose',         false,  {@(v) islogical(v)});
  options = HT_CheckField(options, 'verboseLevel',    0,      {@(v) isscalar(v) && (round(v)==v)});
##  options = HT_CheckField(options, 'swapaxis', false);
  options = HT_CheckField(options, 'checkGeometry',   true,   {@(v) islogical(v)});
  options = HT_CheckField(options, 'strictMatch',     true,   {@(v) islogical(v)});
  options = HT_CheckField(options, 'skipFieldCheck',  false,  @(v) islogical(v));

  if options.verbose, disp(sprintf('Building models <%s> type <%s>', name, lModelTypeName)); tic; endif; % Start timer (tic)

  % Check parameter structure fields. It is desirable to check that no silly errors
  % are made on the parameter name
  if ~options.skipFieldCheck
    assert(all(cellfun(@(v) any(strcmpi(v, {...
      'g', 'T0', 'intersect', 'uv1', 'uv2', 'uv1axis', 'uv2axis', 'epsilon'})),...
       fieldnames(params))), 'Invalid parameter name detected');
  endif

  if iscell(F1)
    assert((numel(F1) == 1) && HT_CheckType(F1, 'face'), 'F1 must be a face structure');
    F1 = F1{1};
  endif

  if iscell(F2)
    assert((numel(F2) == 1) && HT_CheckType(F2, 'face'), 'F2 must be a face structure');
    F2 = F2{1};
  endif

  assert(HT_CheckType(F1, 'face'), 'F1 must be a face structure');
  assert(HT_CheckType(F2, 'face'), 'F2 must be a face structure');

##  assert(isstruct(F1) && isfield(F1, '__type__') && strcmpi(F1.type, 'face'), 'F1 must be a face structure');
  assert(isfield(F1, 'nodes') && iscell(F1.nodes), 'F1.nodes must be a cell array of name');
##  assert(isstruct(F2) && isfield(F2, '__type__') && strcmpi(F2.type, 'face'), 'F2 must be a face structure');
  assert(isfield(F2, 'nodes') && iscell(F2.nodes), 'F2.nodes must be a cell array of name');
  assert(isfield(F1, 'pos') && (size(F1.pos, 2) == 2));
  assert(isfield(F2, 'pos') && (size(F2.pos, 2) == 2));

##  assert(strcmp(class(params.g), 'double') || is_function_handle(params.g), 'Invalid conductance');
  params = HT_CheckField(params, 'g',            Inf,        @(v) isnumeric(v) || is_function_handle(v));
  params = HT_CheckField(params, 'T0',           NaN,        @(v) (numel(v) == 1));
  params = HT_CheckField(params, 'intersect',    false,      {@(v) islogical(v)});
  params = HT_CheckField(params, 'uv1',          [0 1 0 1],  {@(v) (numel(v) == 4) && ~params.intersect});
  params = HT_CheckField(params, 'uv2',          [0 1 0 1],  {@(v) (numel(v) == 4) && ~params.intersect});
  params = HT_CheckField(params, 'uv1axis',      [],         {@(v) (isempty(v) || all(size(v)==[3 2]))});
  params = HT_CheckField(params, 'uv2axis',      [],         {@(v) (isempty(v) || all(size(v)==[3 2]))});
  params = HT_CheckField(params, 'epsilon',      1E-10,       @(v) isnumeric(v) && isscalar(v) && (v>=0));
  assert(~isfield(params, 'C'));

  % Commented 12/09/2023: since the same node could be used by severals polygons of the face
  %assert(numel(unique([F1.nodes(:); F2.nodes(:)])) == numel(F1.nodes)+numel(F2.nodes), 'There is a name conflict between nodes');
  assert(isempty(intersect(F1.nodes, F2.nodes)), 'The faces have nodes in common. It is not supported');
  assert(~xor(isempty(params.uv1axis), isempty(params.uv2axis)));

  uv1 = params.uv1(:)';
  uv2 = params.uv2(:)';

  assert(numel(uv1) == 4 && all(uv1 >= 0) && all(uv1 <= 1));
  assert(numel(uv2) == 4 && all(uv2 >= 0) && all(uv2 <= 1));

  % Check face axis orientation and change them to allow
  % position and distance of nodes to be computed
  if ~isempty(params.uv1axis) && ~isempty(params.uv2axis)
    [F1 uv1] = HT_Face_AdjustAxis(F1, params.uv1axis, uv1, options);
    [F2 uv2] = HT_Face_AdjustAxis(F2, params.uv2axis, uv2, options);
  else
    [F2 uv2] = HT_Face_AdjustAxis(F2, F1, uv2, options);
  endif

  % Depending on the type of faces (based on rectangles only or not)
  % the algorithm differs
  if isempty(F1.vertex) && isempty(F2.vertex)
    M = Int_ConnectFaces_RectRect(M, F1, uv1, F2, uv2, params, options);
  elseif ~isempty(F1.vertex) && ~isempty(F2.vertex)
    assert(columns(F1.dims) == columns(F2.dims), 'Not implemented: F1 and F2 must have the same type of polygons');
    M = Int_ConnectFaces_General(M, F1, uv1, F2, uv2, params, options);
  else
    error('Not implemented: intersection between rect based faces and general faces');
  endif

  M.timer.buildingTime += toc(lTicId);

  if options.verbose, disp(sprintf('Building models <%s> done in %.2f ms', name, 1000*M.timer.buildingTime)); endif; % Start timer (tic)

endfunction

function M = Int_ConnectFaces_General(M, F1, uv1, F2, uv2, params, options)
  HT_ImportConstants();

  assert(options.strictMatch, 'Not implemented. <strictMatch> must be enabled.');
  assert(isscalar(params.g), 'Not implemented. <g> must be a scalar');

  nVertex1 = rows(F1.vertex);
  nVertex2 = rows(F2.vertex);

  % Compute the position in global axis
  lVertex1 = HT_Face_GetAbsoluteData(F1, 'vertex'); %F1.globalPosition + (F1.axis .* F1.size') * F1.vertex';
  lVertex2 = HT_Face_GetAbsoluteData(F2, 'vertex'); %F2.globalPosition + (F2.axis .* F2.size') * F2.vertex';

  % Build a matrix
  lVI = sparse(nVertex1, nVertex2);

  for i=1:nVertex1
##    lVI(i,:) = norm(lVertex2 - lVertex1(:,i), Inf, "cols") < params.epsilon;
    lVI(i,:) = norm(lVertex2 - lVertex1(:,i), Inf, "cols") == 0.0;
  endfor

  nQuads1 = rows(F1.dims);
  nQuads2 = rows(F2.dims);

  lMatchIndex = zeros(nQuads1, 1);

  for i=1:nQuads1
    lSelect = true(nQuads2, 1);
    % subset = Find in F2 all triangles that have the first vertex F1.vertex(i1)        with i1=dims(i,1)
    for k=1:columns(F1.dims)
      lSubSetVertex = int32(find(lVI( F1.dims(i,k), : )));
      if isempty(lSubSetVertex), lSelect = false; break; endif;

      % Find all polygons that use any of the vertex in lSubSetVertex
      t = F2.dims(lSelect,:);
      t = any(t'(:) == lSubSetVertex, 2);
      t = reshape(t, 4, numel(t)/4);
      t = any(t, 1);

      lSelect(lSelect) = lSelect(lSelect) & t';
      if ~any(lSelect), break; endif
    endfor

    lMatch = find(lSelect);
    if isempty(lMatch), continue; endif
    assert(numel(lMatch) == 1, sprintf('Invalid geometry: Multiple polygons (=%d) seem to match', numel(lMatch)));
    lMatchIndex(i) = lMatch;
  endfor

  % Now, all polygons match have been detected
  lMatch1Valid = lMatchIndex != 0;
  lMatch1IndexList = find(lMatch1Valid);
  lMatch1IndexFromFull = zeros(nQuads1, 1);
  lMatch1IndexFromFull(lMatch1IndexList) = 1:sum(lMatch1Valid);

  lMatch2IndexValid = lMatchIndex(lMatch1IndexList);
  lMatch2IndexFromFull = zeros(nQuads2, 1);
  lMatch2IndexFromFull(lMatch2IndexValid) = 1:numel(lMatch2IndexValid);

  lNodesArea = HT_Face_GetAbsoluteData(F1, 'nodesArea');
  lNodesArea(~lMatch1Valid) = [];

  assert(numel(lMatch1IndexList) > 0, sprintf('No intersection found between face <%s> and face <%s>', F1.name, F2.name));

  if options.verbose
    disp(sprintf('Connecting face <%s> and face <%s>. %d matches found. Area=%g', F1.name, F2.name, numel(lMatch1IndexList), sum(lNodesArea)));
  endif

  % Severals occurences of the same node may exist. It is handled below.
  [lNode1Unique, ~, j1] = unique(F1.nodes(lMatch1Valid));
  [lNode2Unique, ~, j2] = unique(F2.nodes(lMatch2IndexValid));

##  M.nodes = [F1.nodes(lMatch1Valid); F2.nodes(lMatch2Valid)];
  M.nodes = [lNode1Unique; lNode2Unique];
  M.G = sparse(numel(lNode1Unique) + numel(lNode2Unique), numel(lNode1Unique) + numel(lNode2Unique));

  lResistance = 1./params.g + F1.r(lMatch1IndexList) + F2.r(lMatch2IndexValid);
  if any(abs(lResistance) < 1E-10)
    M.counter.errors++;
    error(sprintf('Some thermal resistance between <%s> and <%s> is zero', F1.name, F2.name));
  endif

  %find(lMatch1Valid) , lMatchVec(lMatch1Valid), lResistance
  %i1(find(lMatch1Valid)) , numel(lNode1Unique) + i2(lMatchVec(lMatch1Valid)), lResistance
  %  sparse (i, j, sv, 3, 4)

##  i = ind1 ind1 ind2 ind2
##  j = ind2 ind1 ind1 ind2
##  v = +G    -G   +G  -G

##  tic();
  iList = j1(lMatch1IndexFromFull(lMatch1IndexList));
  jList = numel(lNode1Unique)+j2(lMatch2IndexFromFull(lMatch2IndexValid));
  v = lNodesArea ./ lResistance;

  iListTotal = [iList iList jList jList]'(:);
  jListTotal = [jList iList iList jList]'(:);
  v = v .* [1 -1 1 -1];
  v = v'(:);

  nNodes = numel(lNode1Unique) + numel(lNode2Unique);
  M.G = sparse(iListTotal, jListTotal, v, nNodes, nNodes);
  M.T0 = params.T0;
  M.C = NaN(nNodes, 1);

  if isscalar(params.T0), M.T0 = repmat(M.T0, nNodes, 1); endif;

##  disp(sprintf('Time 1 = %f', toc()));

  % Expand version:
  % Start
  % tic();
##  for i=lMatch1IndexList'
##    % Connection is now made between node F1.dims(i) and F2.dims(lMatchPoly)
##    lResistance = 1./params.g + F1.r(i) + F2.r(lMatchIndex(i));
##
##    ind1 = j1(lMatch1IndexFromFull(i));
##    ind2 = numel(lNode1Unique) + j2(lMatch2IndexFromFull(lMatchIndex(i)));
##
##    % Could be probably optimized with sparse matrix construction
##    G = lNodesArea(i) / lResistance;
##    M.G(ind1, ind2 ) += G;
##    M.G(ind1, ind1) -= G;
##    M.G(ind2, ind1) += G;
##    M.G(ind2,ind2) -= G;
##
##    % Check that nodes are correctly connected
####    assert(strcmp(F1.nodes(i), M.nodes(ind1)));
####    assert(strcmp(F2.nodes(lMatchIndex(i)), M.nodes(ind2)));
##  endfor
####  disp(sprintf('Time 2 = %f', toc()));
  % End

endfunction

function M = Int_ConnectFaces_RectRect(M, F1, uv1, F2, uv2, params, options)
  HT_ImportConstants();

  % if option <intersect> is enabled, the intersection of faces is found first
  if params.intersect
    assert(norm(cross(F1.norm, F2.norm), Inf) < 1E-10, "Faces are not oriented the same way. No intersection can be found");
    lDeltaPos = F2.globalPosition - F1.globalPosition;

    lU = [  F1.size(U) * [0 1]; ...
            dot(lDeltaPos, F1.axis(:,U)) + F2.size(U) * [0 1] ];
    lV = [  F1.size(V) * [0 1]; ...
            dot(lDeltaPos, F1.axis(:,V)) + F2.size(V) * [0 1] ];

    lU = [max(lU(:,1)) min(lU(:,2))];
    lV = [max(lV(:,1)) min(lV(:,2))];

    % Now revert the intersection area into the local axis of each faces
    uv1 = [lU ./ F1.size(U), lV ./ F1.size(V)];
    uv2 = [(lU - dot(lDeltaPos, F1.axis(:,U)) ) ./ F2.size(U), (lV - dot(lDeltaPos, F1.axis(:,V)) ) ./ F2.size(V)];

    % Remove numerical errors
    uv1 = floor(uv1/HT_VAR_EPSILON_U + 0.5) * HT_VAR_EPSILON_U;
    uv2 = floor(uv2/HT_VAR_EPSILON_U + 0.5) * HT_VAR_EPSILON_U;

    % Clamp uv1 and uv2 to [0,1]. If there is no intersection, uv1(1) is likely to be greater than uv1(2)
    uv1 = clip(uv1, [0 1]);
    uv2 = clip(uv2, [0 1]);
  endif

  if ((uv1(2)-uv1(1))*(uv1(4)-uv1(3)) < 1E-10) || ((uv2(2)-uv2(1))*(uv2(4)-uv2(3)) < 1E-10)
    warning(sprintf('No intersection between face <%s.%s> and face <%s.%s>', ...
          F1.model, F1.name, F2.model, F2.name));
    return;
  endif

  if isempty(params.uv1axis) && isempty(params.uv2axis) && options.checkGeometry
    F1globalPos = F1.globalPosition + F1.axis * (uv1([1 3])' .* F1.size);
    F2globalPos = F2.globalPosition + F2.axis * (uv2([1 3])' .* F2.size);
    if norm(F1globalPos - F2globalPos) > 1E-10
      M.counter.warnings++;
      warning(sprintf('Invalid geometry between face <%s.%s> at (%f,%f,%f) and face <%s.%s> at (%f,%f,%f)', ...
              F1.model, F1.name, F1globalPos, F2.model, F2.name, F2globalPos));
    endif
    clear F1globalPos F2globalPos;
  endif

  % Check dimensions
  uv1chk = uv1 .* repelem(F1.size', 2); %[F1.size(1) F1.size(1) F1.size(2) F1.size(2)];
  uv2chk = uv2 .* repelem(F2.size', 2); %[F2.size(1) F2.size(1) F2.size(2) F2.size(2)];
  lu1 = (uv1chk(2)-uv1chk(1));
  lv1 = (uv1chk(4)-uv1chk(3));
  lu2 = (uv2chk(2)-uv2chk(1));
  lv2 = (uv2chk(4)-uv2chk(3));
  assert(abs(lu1 - lu2) < 1E-10, sprintf('Wrong faces dimensions. u length = %f and %f', lu1, lu2));
  assert(abs(lv1 - lv2) < 1E-10, sprintf('Wrong faces dimensions. v length = %f and %f', lv1, lv2));
  lThreshold = 1E-10*lu1*lv1; % Remove conductance attached to surface area smaller than this
  clear uv1chk uv2chk lu1 lv1 lu2 lv2;

  % Nodes of face 1 that are concerned by the connection
  nodes1vec = find((F1.dims(:,2) > uv1(1)) & (F1.dims(:,1) < uv1(2)) & (F1.dims(:,4) > uv1(3)) & (F1.dims(:,3) < uv1(4)));
  dims1vec = F1.dims(nodes1vec,:) - [uv1(1) uv1(1) uv1(3) uv1(3)]; % Change origin of coordinates
  dims1vec .*= repelem(F1.size', 2); %[F1.size(1) F1.size(1) F1.size(2) F1.size(2)];
  nodes1resistance = F1.r;
  if isempty(nodes1resistance)
    nodes1resistance = zeros(numel(nodes1vec), 1);
  elseif numel(nodes1resistance) == 1
    nodes1resistance = repmat(nodes1resistance, numel(nodes1vec), 1);
  endif
  dims1vec = floor(dims1vec / HT_VAR_EPSILON_POS + 0.5)*HT_VAR_EPSILON_POS;
  dims1vec = max(dims1vec, 0); % This clamp must be done after nodes1resistance computation

  % Nodes of face 2 that are concerned by the connection
  nodes2vec = find((F2.dims(:,2) > uv2(1)) & (F2.dims(:,1) < uv2(2)) & (F2.dims(:,4) > uv2(3)) & (F2.dims(:,3) < uv2(4)));
  dims2vec = F2.dims(nodes2vec,:) - [uv2(1) uv2(1) uv2(3) uv2(3)];
  dims2vec .*= repelem(F2.size', 2); % [F2.size(1) F2.size(1) F2.size(2) F2.size(2)];
  nodes2resistance = F2.r;
  if isempty(nodes2resistance)
    nodes2resistance = zeros(numel(nodes2vec), 1);
  elseif numel(nodes2resistance) == 1
    nodes2resistance = repmat(nodes2resistance, numel(nodes2vec), 1);
  endif
  dims2vec = floor(dims2vec / HT_VAR_EPSILON_POS + 0.5)*HT_VAR_EPSILON_POS;
  dims2vec = max(dims2vec, 0);

  nNodes = numel(nodes1vec) + numel(nodes2vec);

  % Add nodes to the model
  M.nodes = [F1.nodes(nodes1vec)(:); F2.nodes(nodes2vec)(:)];
  nodes2toMnodes = numel(nodes1vec) + (1:numel(nodes2vec));
  M.G = sparse(nNodes, nNodes);
  M.T0 = params.T0;
  M.C = NaN(nNodes, 1);

  if numel(M.T0) == 1, M.T0 = repmat(M.T0, nNodes, 1); endif;

  if is_function_handle(params.g)
    lgFuncType = HT_GetFunctionType({params.g});
  endif

  % For each nodes in face 1, we look for nodes in face 2 which are in contact
  for i=1:numel(nodes1vec)
    n1ind = nodes1vec(i);
    ldims1 = dims1vec(i,:);  % borders along u,u,v,v axis
    lresistance1 = nodes1resistance(i);

    contactNodes = (dims2vec(:,2) > ldims1(1)) & ...
                   (dims2vec(:,1) < ldims1(2)) & ...
                   (dims2vec(:,4) > ldims1(3)) & ...
                   (dims2vec(:,3) < ldims1(4));
##    nind2vec = find(contactNodes); % Convert logical array to index array

    ldims2 = dims2vec(contactNodes,:);
    if isempty(ldims2), continue; endif

    lresistance2 = nodes2resistance(contactNodes);

    % Surface area between nodes1vec(i) and nodes2
    S1to2vec = [max(ldims2(:,1), ldims1(1)), ...
                min(ldims2(:,2), ldims1(2)), ...
                max(ldims2(:,3), ldims1(3)), ...
                min(ldims2(:,4), ldims1(4))];

    S1to2vec = (S1to2vec(:,2) - S1to2vec(:,1)) .* (S1to2vec(:,4) - S1to2vec(:,3));
    assert(all(S1to2vec >=0));

    lRem = S1to2vec < lThreshold;
    S1to2vec(lRem) = [];
    lresistance2(lRem) = [];
    contactNodes(find(contactNodes)(lRem)) = false;
    clear lRem;

    if isempty(S1to2vec), continue; endif % node1 has no contact with other nodes

    % Add the conductance into the model
    if options.verbose && (options.verboseLevel >= 1) || options.strictMatch
      lStr = ''; %strjoin(F2.nodes{nodes2vec(find(contactNodes))});
      lAreaPart = S1to2vec ./ sum(S1to2vec);

      if options.strictMatch
        assert(all(abs(lAreaPart-1.0) < 1E-10), "Some mismatches were detected between nodes of different faces.");
      endif

      if options.verbose && (options.verboseLevel >= 1)
        n2vec = find(contactNodes);
        for k=1:numel(n2vec) % Iterate over nodes that are in contact with current node of face 1
          lStr = strcat(lStr, sprintf('%s(area=%.2f%% =%.2e m²),', F2.nodes{nodes2vec(n2vec(k))}, 100*lAreaPart(k), S1to2vec(k)));
        endfor
        disp(sprintf('Node %d %s <-> %s %s', nodes1vec(i), F1.nodes{nodes1vec(i)}, ...
                sprintf('%d,', find(contactNodes)), lStr));
        clear n2vec;
      endif
      clear lAreaPart lStr;
    endif


##    node1 = F1.nodes(n1ind);
##    node2vec = F2.nodes(nodes2vec(contactNodes));

    if is_function_handle(params.g)
      % eventually to optimize. The anonymous function stores each time a complete
      % backup of local variables. In that case, it means S1to2vec and lresistance2
      % are put multiple times in memory. Argument needed in the function could be
      % stored in M.Gdyn array and send to the function.
      n2vec = find(contactNodes);
      for k=1:numel(n2vec)
        if lgFuncType == 2
          lgFunc = @(t, i) S1to2vec(k) / (1/params.g(t,i) + lresistance2(k) + lresistance1);
        elseif lgFuncType == 3
          lgFunc = @(T, t, i) S1to2vec(k) / (1/params.g(T, t,i) + lresistance2(k) + lresistance1);
        else
          assert(false, "Invalid function");
        endif

        M.Gdyn = [M.Gdyn; {F1.nodes{n1ind},...
                           F2.nodes{n2vec(k)}, ...
                           lgFunc,...
                           lgFuncType }];
      endfor
    else
      Mind = nodes2toMnodes(contactNodes);

      if abs(1./params.g + lresistance2 + lresistance1) < 1E-10
        lExtStr = '';
        if numel(Mind) > 5, lExtStr = '...'; endif;
        error(sprintf('Thermal resistance between <%s> and <%s%s> is zero', ...
                 M.nodes{i}, strjoin(M.nodes(Mind(1:min(4,end))), ';'), lExtStr));
      endif

      G = S1to2vec ./ (1./params.g + lresistance2 + lresistance1);
      M.G(i, Mind) += G';
      M.G(i, i) -= sum(G);
      M.G(Mind, i) += G;
      for k=1:numel(Mind)
        M.G(Mind(k),Mind(k)) -= G(k);
      endfor
      clear Mind G;
    endif
  endfor

  % Remove empty lines and corresponding nodes
##  lRem = find(max(full(abs(M.G)), [], 2) < lThreshold);
##  if options.verbose, disp(sprintf('Remove nodes %s', sprintf('%d,', lRem))); endif;
##  M.G(lRem, :) = [];
##  M.G(:,lRem) = [];
##  M.T0(lRem) = [];
##  M.nodes(lRem) = [];
##  M.C(lRem) = [];

  % Check conductance matrix
  assert(abs(sum(M.G, 2)) < 1E-10, 'Invalid conductance matrix');
  assert(max(max(abs(M.G - M.G'))) < 1E-10, 'Invalid conductance matrix');

endfunction

