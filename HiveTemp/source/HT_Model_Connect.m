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
% This function defines a thermal connection between two object <obj1> and <obj2>
% It returns a thermal model M with name <name>.
% In the future, this function will be used for all thermal connection. For now
% only the thermal connection between node and face is implemented.
%
% Input parameters:
% .name: [string] the new model name
% .obj1: [object] the first object to connect (node only)
% .obj2: [object] the second object to connect (face only)
% .params: [parameters structure] contains all necessary parameters
%     -> .g: [scalar/column vector/function handle] the conductance (W/K/m²)
%            Specifying a column vector allows to define a space varying conductance. In that
%            case, the number of elements in g must be equal to the number of nodes.
%            Specifying a function handle allows to define a time varying conductance. Moreover
%            it can be space varying as well. It can be defined as @(t,i) or @(T,t,i)
%            "t" is the current time, and "i" the current time index
%            If input argument "T" is used, it means the conductance is non linear
%            "T" will be provided during computations. For a Node to Node connection,
%            "T" will be dimension 2x1 (the temperature for each node will be provided)
%            For nNodes nodes, "T" will be of dimension "nNodes+1" x 1
% .options: [parameters structure] Structure used to specify options.
%     -> .verbose: [boolean] default is false. If true, display information.
%
% Output parameters:
% -> M: [model object] Returns the model associated to this connection
%
function M = HT_Model_Connect(name, obj1, obj2, varargin)
  assert(nargin >= 3, 'Missing input parameters. At least 3 are required');
  assert(ischar(name), 'Invalid input argument <name>. Must be char');
  assert(HT_CheckType(obj1, {'node'}), 'Type of obj1 is invalid');    % Only node is supported for obj1 yet
  assert(all(HT_CheckType(obj2, {'node', 'face'})), 'Type of obj2 is invalid');
  assert(numel(obj1) == 1, 'obj1 must be a single object');

  lTicId = tic();   % Start timer (tic)
  lModelTypeName = 'model_face_connection';

  prop = varargin(1:2:end);
  value = varargin(2:2:end);

  if mod(numel(varargin),2) != 0 % Odd number of parameters ? It means the last
                                 % is a structure which contains options
    assert(isstruct(varargin{end}), 'With a odd number of parameters, last argument of HT_Model_Connect must be structure');
    lOptions = varargin{end};
    lOptions = CheckField(lOptions, 'verbose', false, @(v) (numel(v)==1) && islogical(v));
  else
    lOptions = struct(  'verbose', false);
  endif

  lParams = struct('obj1', obj1, 'obj2', obj2); % Store input arguments into structure params used to keep a track of the history
  M = HT_Model_Init(lModelTypeName, name, lParams);
  clear lParams;

  if lOptions.verbose, disp(sprintf('Building models <%s> type <%s>', name, lModelTypeName)); endif;

  % This function is capable to handling a cell array of obj2
  % but obj1 must be a single object
  if iscell(obj2)
    % Make sure all obj2 have the same type
    assert(numel(unique(cellfun(@(v) v.__type__, obj2, 'UniformOutput', false))) == 1, 'All obj2 objects must have the same type');
  else
    % Make sure all obj2 have the same type
    assert(numel(unique(arrayfun(@(v) v.__type__, obj2, 'UniformOutput', false))) == 1, 'All obj2 objects must have the same type');
  endif


  if HT_CheckType(obj1, 'node')
    M = Int_ConnectNode(M, obj1, obj2, prop, value, lOptions);
  elseif HT_CheckType(obj1, 'face')
    error('Not implemented');
  endif

  % Check conductance matrix
  assert(abs(sum(M.G, 2)) < 1E-10, 'Invalid conductance matrix');
  assert(max(max(abs(M.G - M.G'))) < 1E-10, 'Invalid conductance matrix');

  M.timer.buildingTime += toc(lTicId);

  if lOptions.verbose, disp(sprintf('Building models <%s> done in %.2f ms', name, 1000*M.timer.buildingTime)); endif; % Start timer (tic)

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
      assert(chkFunc{i}(getfield(P, field)))
    endif
  endfor

endfunction

% Check the type of a function handle
% time varying or non linear time varying
function M = Int_ConnectNode(M, obj1, objList, prop, value, lOptions)
  if lOptions.verbose, disp(sprintf('Connecting node <%s> to %d object(s)', obj1.name, numel(objList))); endif; % Start timer (tic)

  nObj = numel(objList);
  lConductance = [];

  for i=1:numel(prop)
    if strcmpi(prop{i}, 'g')
      if isnumeric(value{i})
        lConductance = { value{i} };
      elseif is_function_handle(value{i})
        lConductance = { value{i} };
      elseif iscell(value{i})
        assert(numel(value{i}) == nObj, 'Invalid parameter <g>. Cell size is not valid.');
        lConductance = value{i};
      else
        error('Field <g> has not a valid value');
      endif

      if numel(lConductance) == 1, lConductance = repmat(lConductance, nObj,1); endif;
      % Check values
      assert(all(cellfun(@(v) isnumeric(v) || (HT_GetFunctionType(v) > 0), lConductance)), 'Field <g> has not a valid value');
    endif
  endfor

  assert(~isempty(lConductance), 'Field conductance <g> is not defined');
  assert(strcmpi(obj1.mode, 'distributed'), sprintf('Node <%s>.mode must be distributed', obj1.name));

  M.nodes = {obj1.name};
  M.G = sparse(1,1);
  M.T0 = NaN;
  M.C = NaN;

  % Connection between 1 node and somes faces
  if all(HT_CheckType(objList, 'face'))
    for i=1:nObj
      if isstruct(objList)
        obj = objList(i);
      else
        obj = objList{i};
      endif
      lNodesIndex = []; % Store the position of each nodes in this face in the thermal matrices M.G, M.C
      nNodes = numel(M.nodes);

      % Retrieve the list of nodes and add them to the node list
      lNewNodes = setdiff(obj.nodes, M.nodes);    % Remove nodes that are already present in M.nodes
      M.nodes = [M.nodes; lNewNodes];    % lNewNodes = obj.nodes(ia)

##      if numel(lNewNodes) == numel(obj.nodes)  % All nodes were new ?
##        s_idx = nNodes + ia; %((nNodes+1):(nNodes+numel(lNewNodes)))';
##        tf = true(size(s_idx));
##      else % Otherwise, it is slower since we have to find in the model node list the index of each nodes
        % Overwrite <ia> in this case. It indicates to which index in obj.nodes corresponds the index in M.nodes
##        [~, ia, lNodesIndex] = intersect(obj.nodes, M.nodes);
        % 13/09/2023: replaced intersect by ismember since obj.nodes may contain several times the same node
        [tf, s_idx] = ismember(obj.nodes, M.nodes);  % obj.nodes(tf) = M.nodes(s_idx(tf))
##      endif

      nNodes = numel(M.nodes);

      M.G = resize(M.G, nNodes, nNodes);  % Increase the size of the matrix to include new nodes
      M.T0 = postpad(M.T0, nNodes, NaN);
      M.C = postpad(M.C, nNodes, NaN);

      % Fill the matrix with the appropriate coefficient
      if isnumeric(lConductance{i}) && (size(lConductance{i}, 2) == 1) % Column vector ?
        % Retrieve the resistance density [Km²/W] and the area [m²] for each node of the face
        % Sum the resistance of the face + the conductance specified in this thermal connection
        lResistance = (HT_Face_GetAbsoluteData(obj, 'r') + 1./lConductance{i}) ./ HT_Face_GetAbsoluteData(obj, 'nodesArea') ...
                        + obj1.resistance; % Add the absolute resistance defined by the node (generally 0)
        % The nodes order returned by setdiff above may not be identical to obj.nodes
        % so nodes must be rearranged as specified in <ia>
##        lG = 1./lResistance(ia);
        lG = 1./lResistance(tf);
        iList = s_idx(tf);
        jList = ones(size(lG));

        iListTotal = [iList iList jList jList]'(:);
        jListTotal = [jList iList iList jList]'(:);
        lG = lG .* [1 -1 1 -1];
        lG = lG'(:);

        % obj.nodes may contain several times the same node
        % The function sparse implicitely sum all nodes
        M.G += sparse(iListTotal, jListTotal, lG, nNodes, nNodes);

##  iList = j1(lMatch1IndexFromFull(lMatch1IndexList));
##  jList = numel(lNode1Unique)+j2(lMatch2IndexFromFull(lMatch2IndexValid));
##  v = lNodesArea ./ lResistance;
##
##  iListTotal = [iList iList jList jList]'(:);
##  jListTotal = [jList iList iList jList]'(:);
##  v = v .* [1 -1 1 -1];
##  v = v'(:);

##        M.G(1, lNodesIndex) += lG';
##        M.G(lNodesIndex, 1) += lG;
##        lDiag = diag(M.G);
##        lDiag(1) -= sum(lG);
##        lDiag(lNodesIndex) -= lG;
##        M.G = spdiags(lDiag,0, M.G); % Replace the diagonal coefficients
##        clear lG lDiag lResistance;
      elseif isnumeric(lConductance{i}) && (size(lConductance{i}, 2) > 1) % Varying with time ? (severals columns)
        assert(any(size(lConductance{i},1) == [1 nNodes]), 'Matrix for non uniform conductance must have 1 or <nNodes> lines');

        lArea = HT_Face_GetAbsoluteData(obj, 'nodesArea');
        lResistance = (1./lArea) .* (HT_Face_GetAbsoluteData(obj, 'r') + 1./lConductance{i})  ...
                        + obj1.resistance; % Add the absolute resistance defined by the node (generally 0)
        M.Gdyn = [M.Gdyn; { obj1.name, ...
                            obj.nodes, ...
                            1./lResistance, ...
                            0}];
      elseif is_function_handle(lConductance{i})
        lgFuncType = HT_GetFunctionType(lConductance{i});
        lgFunc = lConductance{i};
        lArea = HT_Face_GetAbsoluteData(obj, 'nodesArea');
        lRDensity = HT_Face_GetAbsoluteData(obj, 'r');

        if lgFuncType == 2
          lR = obj1.resistance;

          M.Gdyn = [M.Gdyn; { obj1.name, ...
                              obj.nodes, ...
                              @(t, i) 1 ./ ( (lRDensity + 1/lgFunc(t,i)) ./ lArea + lR), ...
                              lgFuncType}];
          clear lR;
        elseif lgFuncType == 3
          lR = obj1.resistance;

          M.Gdyn = [M.Gdyn; { obj1.name, ...
                              obj.nodes(k), ...
                              @(T, t, i) 1 ./ ( (lRDensity + 1/lgFunc(T, t,i)) ./ lArea + lR), ...
                              lgFuncType}];
          clear lR;
        else
          assert(false, "Invalid function");
        endif

      else
        error('Invalid value of conductance');
      endif
    endfor

  else
    error('Invalid type of objList');
  endif
endfunction
