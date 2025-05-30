%  This file is part of project HiveTemp.
%  This work was supported by the Better-B project, which has received funding
%  from the European Union, the Swiss State Secretariat for Education, Research
%  and Innovation (SERI) and UK Research and Innovation (UKRI) under the UK
%  government's Horizon Europe funding guarantee (grant number 10068544). Views
%  and opinions expressed are however those of the author(s) only and do not
%  necessarily reflect those of the European Union, European Research Executive
%  Agency (REA), SERI or UKRI. Neither the European Union nor the granting
%  authorities can be held responsible for them.
%
%  Copyright (c) 2025 AltRD-Emmanuel Ruffio
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
% Résolution du modèle
% Chaque élément de commande (liste de structure <C>)
% .name = nom de la commande
% .data = scalar/matrix/function handle @(t,i,[index],[userData])
%         If <index> is present, it specifies a non uniform type command.
%         If an array of face structure is specified as .nodes field, the attached
%         function handle will be called for each structure independently.
%         Function handle must then have the following definition
%         function handle @(t,i,index,face,[userData])
%         For non uniform flux:
%         if function handle returns a cell array, it must have 2 elements:
%         1) The first is a scalar will go in U vector
%         2) The second is a vector that will go in command matrix B
%
% Arguments d'entrée
% .T0 = température initiale [1x1] ou [1xnNode]
% .tVec = vecteur temps
%         ou cell dim3x1 = {startTime, dt, nt}
% .options = structure d'options
%     -> .verbose
%     -> .info : 'none', 'iter', ['default'] Print some information during computations
%     -> .all : return all node temperature at all times  T = dim nNode x nTimes
%     -> .replaceT0 : does not take into account previous specification about initial temperature
%                     and set all initial temperature to T0
%     -> .algorithm : 'exact'
%     -> .mergeTU : [true] add nodes whose temperature is commanded into the T matrix and Gnodes array
%     -> .maxdt : [scalar] maximum value of the time step. If the specified time vector
%                 contains too large time step, it is divided into small times steps
%     -> .iterstartfunc : function handle @(k,t,Told, Uold) a user defined
%                         function which returns <userData> structure that will
%                         be passed to command function <cmd.data> at the last
%                         argument. <iterstartfunc> can be used to do some
%                         pre computation that is used several times by other
%                         cmd.data function (ex: sun position)
%                         FOr the first iteration, Uold is empty.
%     -> .unit = '', 'degres', 'kelvin': specify the unit used to specify temperature.
%                         If radiation problems are specified, the proper conversion will be done.
%
% arguments de sortie
% T vecteur des températures [(1 ou nTimes) x nNode]
% U vecteur des commandes [(1 ou nTimes) x nCmd]
% Gnodes le nom des noeuds de chaque ligne du vecteur T

function [Tmat, Umat, Gnodes, Unodes, t] = HT_SolveModel(M, Cmd, T0, tVec, options)
  HT_ImportConstants();

  assert(nargin >= 4, 'There must be at least 4 input arguments');
  lRunningTime = cputime();
  lRunningTimeId = tic();

  if options.verbose, disp(sprintf('Solving model <%s>', M.name)); tic; endif; % Start timer (tic)

  if nargin < 5, options = struct(); endif;
  assert(all(cellfun(@(v) sum(strcmpi(v, ...
    {'verbose', 'info', 'all', 'replaceT0', 'algorithm', 'mergeTU', 'iterstartfunc', 'maxdt', 'unit'})), fieldnames(options)) == 1), 'Invalid options structure');
  options = HT_CheckField(options, 'verbose', 4);
  options = HT_CheckField(options, 'info', 'default', {@(v) any(strcmpi(v, {'none', 'iter', 'default'}))});
  options = HT_CheckField(options, 'all', false);
  options = HT_CheckField(options, 'replaceT0', false);
  options = HT_CheckField(options, 'algorithm', 'exact');
  options = HT_CheckField(options, 'mergeTU', true);
  options = HT_CheckField(options, 'maxdt', []);
  options = HT_CheckField(options, 'iterstartfunc', []);
  options = HT_CheckField(options, 'unit', [], @(v) any(strcmpi(v, {'', 'degres', 'kelvin'})));

  assert(any(strcmpi(options.algorithm, {'exact', 'linear', 'lu', 'chol'})), sprintf('Invalid algorithm <%s>', options.algorithm));
  assert(isempty(options.iterstartfunc) || (numel(Int_GetFunctionArguments(options.iterstartfunc)) == 4), ...
    "<.iterstartfunc> must have 4 arguments: k,t,Told, Uold");

  if options.verbose
    if strcmpi(options.algorithm, 'exact')
      disp('Using exact method: minimisation of residuals');
    elseif strcmpi(options.algorithm, 'linear')
      disp('Using linearization method');
    elseif strcmpi(options.algorithm, 'lu')
      disp('Using LU decomposition');
    elseif strcmpi(options.algorithm, 'chol')
      disp('Using cholestky decomposition');
    endif
  endif

  assert((iscell(tVec) && (numel(tVec)==3)) || isnumeric(tVec), 'Invalid time vector');
  if ~iscell(tVec)
    tStart = tVec(1);
    tEnd = tVec(end);
    nt = numel(tVec);

    dt = tVec(2:end)-tVec(1:(end-1));
    mindt = min(dt);
    maxdt = max(abs(dt - mindt)/mindt);
    if maxdt < 1E-10
      warning("HT:Analysis", sprintf('Very small variations of time step (=%.1E) was found. Consider using t={startTime, finalTime, dt} specification to improve computation speed', maxdt));
    endif
    clear dt mindt maxdt;
  elseif iscell(tVec)
    tStart = tVec{1};
    tEnd = tVec{2} * tVec{3};
    nt = tVec{3};
  endif

  % Radiation problems is present ?
  if ~isempty(M.rad)
    if isempty(options.unit) && (any(M.T0 < 250) || any(T0 < 250))
      warning("HT:Analysis", 'Temperature below 200K and some radiation problems were detected. Please check the temperature unit.');
    endif
  endif

  %
  lConvertToKelvin = ~isempty(M.rad) && strcmpi(options.unit, 'degres');
  if lConvertToKelvin
    lConvertUnitOffset = 273.15;
  else
    lConvertUnitOffset = 0;
  endif

  % Vérification des commandes
  assert(all(HT_CheckType(Cmd, 'command')), 'Invalid command list');
  lCmdNameList = cellfun(@(v) v.name, Cmd, "UniformOutput", false);
  % Duplicate command ?
  assert(numel(unique(lCmdNameList)) == numel(lCmdNameList), 'A node in command vector must appear only once');
  clear lCmdNameList;

  % And specify their type
  % .sys=11 : temperature constant and uniform
  % .sys=12 : temperature variable and uniform
  % .sys=13 : interpolated temperature variable and uniform
  % .sys=14 : (function(k) temperature variable and uniform
  % .sys=15 : (function(t) temperature variable and uniform
  % .sys=21 : flux constant and uniform
  % .sys=22 : variable and uniform flux
  % .sys=23 : constant and non uniform flux
  % .sys=24 : variable and non uniform flux
  % .sys=25 : interpolated variable and non uniform flux
  % .sys=26 : (function)(k) variable and uniform flux
  % .sys=27 : (function)(t) variable and uniform flux
  % .sys=28 : (function)(t,k,index) variable and uniform flux

  if options.verbose >= 4, disp(sprintf('  -> Analyzing command list')); endif;


  % ----- First pass ------ is devoted to temperature type command
  % This kind of command removes nodes from the model since their temperature does
  % not have to be solved anymore. We keep track of all removed nodes.
  nCmd = numel(Cmd);
  lNodeTempCmdList = cell(nCmd, 1);     % Store for each command the list of node to remove
  B = sparse(numel(M.nodes), 1);        % Flux vector is prepared as well

  % Construct a vector containing the index of all temperature type command
  lTemperatureCmdIndex = cellfun(@(v) strcmpi(v.type, 'temperature'), Cmd);
  lTemperatureCmdIndex = uint32(find(lTemperatureCmdIndex));

  for i=int32(1:nCmd)
    cmd = Cmd{i};

    % Append internal field to command structure
    cmd.sys = struct( 'type', 0,                ... % Command type: flux/temperature/uniform flux/constant flux...
                      'nodeIndex', [],          ... % List of valid nodes
                      'nodeFlag', [],
                      'nodeRemoved', [],        ... %
                      'BT', [],                 ... % Command matrix for temperature type command
                      'hasUserData', false,     ... % Specify if <userData> struct must be sent to the function
                      'faceProvided', false,    ... % Specify if the user gave list of faces or list of nodes
                      'faceIndex', {{}},        ... % Cell [dim nf x 1] that contains integer index of each valid nodes
                      'faceIndexUnique', {{}},  ... % Similar to .faceIndex but contains only one occurence of each nodes.
                      'faceIndexOperation', {{}},   ... % It is the result of [~, ~, faceIndexOperation] = unique(faceIndex). It specifiy in will column of BT matrix each nodes is accumulated
                      'faceFlag', {{}},         ... % Face flag does not have the same size as faceIndex since it contains the old nodes (before removal operation)
                      'faceRemoved', []         ... % Contains the index in Uvec (removed node with temperature command) of every index in command (faceRemoved == 0) for index that were not removed
                      );
    lInfos = HT_Cmd_GetType(cmd, options);
    cmd.sys.type = lInfos.type;
    cmd.sys.hasUserData = lInfos.userData;

    if options.verbose, disp(lInfos.typeString); endif;

    if any(HT_CheckType(cmd.nodes, 'face'))
      assert(all(HT_CheckType(cmd.nodes, 'face')));
      nFace = numel(cmd.nodes);

      cmd.sys.faceIndex = cell(nFace, 1);
      cmd.sys.faceRemoved = cell(nFace, 1); % Only used with temperature type command
      cmd.sys.faceFlag = cell(nFace, 1);
      for k=1:nFace
        [~, s_idx] = ismember(cmd.nodes{k}.nodes, M.nodes);
        assert(all(s_idx > 0), sprintf('%d invalid command node(s) found <%s>. It does not exist in the model', numel(s_idx == 0), strjoin(cmd.nodes{k}.nodes(s_idx == 0), ',')));
        cmd.sys.faceIndex{k} = uint32(s_idx);
      endfor

      cmd.sys.faceProvided = true;
      cmd.sys.BT = cell(nFace, 1);
    else
      [~, s_idx] = ismember(cmd.nodes, M.nodes);
      assert(s_idx > 0, sprintf('Invalid command node found <%s>. It does not exist in the model', cmd.nodes{1}));
      cmd.sys.nodeIndex = uint32(s_idx);
    endif

    % ========== Temperature type ===========
    if strcmpi(cmd.type, 'temperature')

      if cmd.sys.faceProvided

        for k=1:numel(cmd.nodes) % Loop over faces
          [lFaceIndexUnique, ~, lFaceIndex_j] = unique(cmd.sys.faceIndex{k});
          cmd.sys.BT{k} = M.G(:, lFaceIndexUnique );
          % Save the merge index list to be used during each time step
          cmd.sys.faceIndexUnique{k} = lFaceIndexUnique;
          cmd.sys.faceIndexOperation{k} = lFaceIndex_j;

          if lConvertToKelvin, cmd.data{k} += lConvertUnitOffset; endif

          % If temperature is constant (uniform or not uniform), it is applied now to the command matrix B
          if cmd.sys.type(k) == HT_BTYPE_T_CONSTANT_UNIFORM
            B += sum(cmd.sys.BT{k},2) * cmd.data{k};
          elseif cmd.sys.type(k) == HT_BTYPE_T_CONSTANT_NON_UNIFORM
##            B += cmd.sys.BT{k} * cmd.data{k}; % <= Old version incompatible with circular face since severals identical nodes are specified in the face
            % Node areas are necessary here when nodes may be specified several times in a face
            % In that case, the temperature on duplicate nodes is computed using
            % weighted average algorithm. Function accumdim is used to accumulate
            % the contribution of each occurence of duplicate nodes.
            t = accumdim(lFaceIndex_j, [cmd.data{k} .* cmd.area{k} , cmd.area{k}]);
            t = t(:,1) ./ t(:,2);
            B += cmd.sys.BT{k} * t;
          elseif cmd.sys.type(k) == HT_BTYPE_T_VARIABLE_UNIFORM
            % If temperature is uniform, we can sum all columns of BT into a single one
            % But it not inserted into command vector B since it is variable.
            cmd.sys.BT{k} = sum(cmd.sys.BT{k}, 2);
          elseif cmd.sys.type(k) == HT_BTYPE_T_VARIABLE_NON_UNIFORM
            % Like CONSTANT_NON_UNIFORM, the temperature on duplicate nodes is computed using
            % weighted average algorithm.
            t = accumdim(lFaceIndex_j, [cmd.data{k} .* cmd.area{k} , cmd.area{k}]);
            t = t(:,1) ./ t(:,2);

          elseif cmd.sys.type(k) == HT_BTYPE_T_INTERPOLATION_ARRAY % Two columns data: [Temperature, Time]
            % In this case, the temperature is variable and uniform
            % so each columns can be merged.
            cmd.sys.BT{k} = sum(cmd.sys.BT{k}, 2);
          else
            error(sprintf('SolveModel::Temperature boundary type <%d> is Not implemented', cmd.sys.type(k)));
          endif
        endfor

        % A list of node index is made for each temperature command type
        % It is legal to have the same node referenced several times in the same
        % command, but it is illegal if the same node is referenced by different
        % commands. To avoid an error during the check of unicity in lNodeTempCmdList
        % below, the set of nodes using <unique> function is done here.
        ##        lNodeTempCmdList{i} = uint32(cell2mat(cmd.sys.faceIndex));
        lNodeTempCmdList{i} = unique(cell2mat(cmd.sys.faceIndex));
      else
        lNodeTempCmdList{i} = uint32(cmd.sys.nodeIndex);

        cmd.sys.BT = M.G(:, cmd.sys.nodeIndex );

        if lConvertToKelvin, cmd.data += lConvertUnitOffset; endif

        % If the flux is constant (uniform or not uniform), it is applied now to the command matrix B
        if cmd.sys.type == HT_BTYPE_T_CONSTANT_UNIFORM
          B += sum(cmd.sys.BT, 2) * cmd.data;
        elseif cmd.sys.type == HT_BTYPE_T_CONSTANT_NON_UNIFORM
          B += cmd.sys.BT * cmd.data;
        elseif cmd.sys.type == HT_BTYPE_T_VARIABLE_UNIFORM
          % If temperature is uniform, we can sum all columns of BT into a single one
          cmd.sys.BT = sum(cmd.sys.BT, 2);
        elseif cmd.sys.type(k) == HT_BTYPE_T_VARIABLE_NON_UNIFORM
          % Nothing to do
        elseif cmd.sys.type(k) == HT_BTYPE_T_INTERPOLATION_ARRAY % Two columns data: [Temperature, Time]
          % In this case, the temperature is variable and uniform
          % so each columns can be merged.
          cmd.sys.BT = sum(cmd.sys.BT, 2);
        else
          error(sprintf('SolveModel::Temperature boundary type <%d> is Not implemented', cmd.sys.type(k)));
        endif
      endif
    else % ========== Flux type ===========

      if cmd.sys.faceProvided
        % If multiple faces are specified, some nodes may be concerned by severals faces
        % The contribution of each faces must be added
        for k=1:numel(cmd.sys.faceIndex)
          % If the flux is constant (uniform or not uniform), it is applied now to the command matrix B
          if any(cmd.sys.type(k) == [HT_BTYPE_FLUX_CONSTANT_UNIFORM   HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM])
            if isempty(cmd.area) || isempty(cmd.area{k})
              B( cmd.sys.faceIndex{k}) += cmd.data{k};
            else
##              B( cmd.sys.faceIndex{k}) += cmd.data{k} .* cmd.area{k};
              B += accumarray(cmd.sys.faceIndex{k}, cmd.data{k} .* cmd.area{k}, size(B));
            endif
          endif
        endfor

      else

        % If the flux is constant (uniform or not uniform), it is applied now to the command matrix B
        if any(cmd.sys.type == [HT_BTYPE_FLUX_CONSTANT_UNIFORM   HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM])
          if isempty(cmd.area)
            B( cmd.sys.nodeIndex) += cmd.data;
          else
            B( cmd.sys.nodeIndex) += cmd.data .* cmd.area;
          endif
        endif
      endif
    endif

    Cmd{i} = cmd;
  endfor

  clear cmd lInfos;

  lNodeTempCmdList = cell2mat(lNodeTempCmdList);
  nNodesRemoved = numel(lNodeTempCmdList);
  % Make sure that each node is commanded once among temperature type command
  assert(numel(unique(lNodeTempCmdList)) == numel(lNodeTempCmdList), 'Some nodes have multiple temperature command');

  % Display a warning if flux type command are applied to nodes with temperature type commmand
  % it is not necessarily an error
  if ~all(~cellfun(@(v) strcmpi(v.type, 'flux') && any(ismember(v.sys.nodeIndex, lNodeTempCmdList)), Cmd))
    warning('HT:Inconcistency', 'Some nodes with flux type boundary condition are also subjected to temperature type boundary condition');
  endif

  lNodeRemoveFlag = false(numel(M.nodes), 1);
  lNodeRemoveFlag(lNodeTempCmdList) = true;

  t = (1:numel(M.nodes))';
  t(lNodeTempCmdList) = [];
  [tf s_idx] = ismember((1:numel(M.nodes))', t);
  lNodeIndexUpdate = NA(numel(M.nodes), 1);
  lNodeIndexUpdate(tf) = s_idx(tf);

  nNodes = numel(M.nodes) - numel(lNodeTempCmdList);

  % Allocate temperature matrix
  if options.all
    Tmat = zeros(nNodes, nt);
  else
    Tmat = zeros(nNodes, 1);
  endif

  % If option mergeTU is true, an other matrix will store temperature command value
  if options.mergeTU
    Umat = zeros(nNodesRemoved, columns(Tmat));
    Unow = zeros(nNodesRemoved, 1);
    Unodes = M.nodes(lNodeTempCmdList);
  else
    Umat = [];
    Unow = [];
  endif

  % 2nd pass: every nodes that has to be removed is removed.
  for i=int32(1:nCmd)
    if strcmpi(Cmd{i}.type, 'temperature')

      if Cmd{i}.sys.faceProvided
        for k=1:numel(Cmd{i}.nodes)
          % Store in faceRemoved
          % First output argument (tf) is not necessary here since all nodes of temperature
          % type command were removed. In other words, tf is true for every index.
          [~, Cmd{i}.sys.faceRemoved{k}] = ismember(Cmd{i}.sys.faceIndexUnique{k}, lNodeTempCmdList);
          Cmd{i}.sys.BT{k}(lNodeTempCmdList,:) = [];
          Cmd{i}.sys.faceIndex{k} = lNodeIndexUpdate(Cmd{i}.sys.faceIndex{k});
##          Cmd{i}.sys.faceIndexUnique{k} = lNodeIndexUpdate(Cmd{i}.sys.faceIndexUnique{k});

          if options.mergeTU
            if Cmd{i}.sys.type == HT_BTYPE_T_CONSTANT_UNIFORM
              Umat( Cmd{i}.sys.faceRemoved{k} ,:) = Cmd{i}.data{k};
            elseif Cmd{i}.sys.type == HT_BTYPE_T_CONSTANT_NON_UNIFORM
##              Umat(Cmd{i}.sys.faceRemoved{k},:) = Cmd{i}.data{k};

              t = accumdim(Cmd{i}.sys.faceIndexOperation{k}, [Cmd{i}.data{k} .* Cmd{i}.area{k} , Cmd{i}.area{k}]);
              Umat(Cmd{i}.sys.faceRemoved{k},:) = repmat(t(:,1) ./ t(:,2), 1, columns(Umat));
            endif
          endif

        endfor
      else
        [~, Cmd{i}.sys.nodeRemoved] = ismember(Cmd{i}.sys.nodeIndex, lNodeTempCmdList);
        Cmd{i}.sys.BT(lNodeTempCmdList,:) = [];

        if options.mergeTU
          if Cmd{i}.sys.type == HT_BTYPE_T_CONSTANT_UNIFORM
            Umat(Cmd{i}.sys.nodeRemoved,:) = Cmd{i}.data;
          elseif Cmd{i}.sys.type == HT_BTYPE_T_CONSTANT_NON_UNIFORM
            Umat(Cmd{i}.sys.nodeRemoved,:) = Cmd{i}.data;
          endif
        endif
      endif

    else % If strcmpi(Cmd{i}.type, 'temperature')
      if Cmd{i}.sys.faceProvided
        for k=1:numel(Cmd{i}.nodes)
          lRemoveList = lNodeRemoveFlag(Cmd{i}.sys.faceIndex{k});
          Cmd{i}.sys.faceIndex{k}(lRemoveList) = [];
          Cmd{i}.sys.faceFlag{k} = lRemoveList;
          Cmd{i}.sys.faceIndex{k} = lNodeIndexUpdate(Cmd{i}.sys.faceIndex{k});

          if rows(Cmd{i}.data{k}) > 1
            Cmd{i}.data{k}(lRemoveList, :) = [];
          endif

          if rows(Cmd{i}.area{k}) > 1
            Cmd{i}.area{k}(lRemoveList, :) = [];
          endif

          if isempty(Cmd{i}.sys.faceIndex{k})
            warning('HT:Inconsistency', sprintf('Cmd <%s.face(%s)> has become ineffective. All nodes were removed', Cmd{i}.name, Cmd{i}.nodes{k}.name));
          endif
        endfor
      else
        lRemoveList = lNodeRemoveFlag(Cmd{i}.sys.nodeIndex);
        Cmd{i}.sys.nodeIndex(lRemoveList) = [];
        Cmd{i}.sys.nodeFlag = lRemoveList;
        Cmd{i}.sys.nodeIndex = lNodeIndexUpdate(Cmd{i}.sys.nodeIndex);

        % Remove all rows of data that is related to a removed node
        if rows(Cmd{i}.data) > 1
          Cmd{i}.data(lRemoveList,:) = [];
        endif

        if rows(Cmd{i}.area{k}) > 1
          Cmd{i}.area{k}(lRemoveList, :) = [];
        endif

        if isempty(Cmd{i}.sys.nodeIndex)
          warning('HT:Inconsistency', sprintf('Cmd <%s> has become ineffective. All nodes were removed', Cmd{i}.name));
        endif
      endif
    endif % If strcmpi(Cmd{i}.type, 'temperature')

  endfor

  % Apply the remove operation to model matrices
  lKeepNodeFlag = ~lNodeRemoveFlag;

  Gnodes = M.nodes(lKeepNodeFlag);

  G = M.G(lKeepNodeFlag,lKeepNodeFlag);   % Remove the equation corresponding to each temperature command
  B(lNodeTempCmdList) = [];
  C = M.C(lKeepNodeFlag);                 % Same operation for capacity matrix
  assert(~any(isnan(C)), sprintf('Some capacitance are not defined (%d nodes: %s)', sum(isnan(C)), HT_StrJoin(Gnodes(isnan(C)), '; ', 5)));

  C = spdiags(C, 0, nNodes, nNodes); % Convert vector to diagonal capacity matrix


  % Selection matrix and its transpose, and R matrix = A(F-Id)(Id-(Id-epsilon)F) epsilon sigma
  lRadiationMatrix = INT_PrepareRadiationMatrix(M.rad, M.nodes); % Returns the selection matrix Sk for each radiation problem
  lRadiationMatrixCmd = lRadiationMatrix(lKeepNodeFlag, lNodeTempCmdList);
  lRadiationMatrix = lRadiationMatrix(lKeepNodeFlag,lKeepNodeFlag);

  % Matrix lRadiationMatrixCmd is empty ? It means that none of the nodes that are
  % part of a radiation problem are in the temperature command vector
  % In that case, it is cleared.
  ind = any(lRadiationMatrixCmd);   % Look for non empty columns
  if any(ind)
    lNodeTempRadiationList = ind;
    lRadiationMatrixCmd = lRadiationMatrixCmd(:,ind);
  else
    lNodeTempRadiationList = [];
    lRadiationMatrixCmd = [];
  endif

  % If options.replaceT0 is true, all intial temperatures are replaced by the specified vector T0
  assert(isscalar(T0) || (isempty(T0) && options.replaceT0), 'Invalid initial temperature. Must be a scalar value.');
  if ~options.replaceT0
    % In that case, only NaN values of the model initial temperature are set to T0
    tmp = M.T0(lKeepNodeFlag);
    ind = isnan(tmp);
    assert(isempty(ind) || ~isempty(T0), ...
      sprintf('Some initial temperature are not defined and T0 is empty. Following nodes:\n%s', HT_StrJoin(Gnodes(ind), ';', 10)));
    tmp(ind) = T0;
    T0 = tmp;
    clear tmp ind;
  endif

  assert(issymmetric(G), 'Conductance matrix is not diagonal');
  tmp = sum(G,2);

  for i=1:numel(lTemperatureCmdIndex)
    ind = lTemperatureCmdIndex(i);

    if strcmpi(Cmd{ind}.type, 'temperature')
      if Cmd{ind}.sys.faceProvided
        for k=1:numel(Cmd{ind}.nodes)
          tmp += sum(Cmd{ind}.sys.BT{k},2);
        endfor
      else
        tmp += sum(Cmd{ind}.sys.BT, 2);
      endif
    endif
  endfor
##  tmp = abs(sum(G,2)); % + sum(B(:, cellfun(@(v) strcmpi(v.type, 'temperature'), Cmd) ),2));
  tmp = abs(tmp);
  if options.verbose, disp(sprintf('Max Conductance matrix residuals is %E', max(tmp))); endif;
  tmp = (abs(tmp) ./ abs(diag(G))) < 1E-12;
  assert(all(tmp), sprintf('Some line (%s) are invalid: sum != 0', HT_StrJoin(int32(find(full(~tmp), 10)), ',', 10)));

  % Time varying conductance
  % To speed up update operation, M.Gdyn cell arrays is changed by
  % replacing nodes names in column 1 and 2 by node ids
  % M.Gdyn{i,1:4} contains (nodes1, nodes2, function, function type)
  % nodes2 can be a cell array of nodes
  % Gdyn{i,1:5} will contains (nodes1, nodes2, cmdId, function, function type)
  %             cmdId is the index in command array
  nDyn = size(M.Gdyn,1);
  Gdyn = cell(nDyn, 5); % Columns 1,2 for conductance matrix, and 3 for command matrix column, and function type

  if ~isempty(M.Gdyn)
    lNodesToValidNodes = 1:numel(M.nodes);
    lNodesToValidNodes(lNodeTempCmdList) = 0;
    lNodesToTCmdNodes = zeros(numel(M.nodes), 1);
    lNodesToTCmdNodes(lNodeTempCmdList) = 1:numel(lNodeTempCmdList);

    for i=1:nDyn
      if M.Gdyn{i,4} == 0 % Function type == 0 means a matrix is provided instead of a function
        assert(any(columns(M.Gdyn{i,3}) == [1 nt]), 'Matrix for time varying conductance must have 1 or <nt> columns');
      endif

      % Find the node in node list
      % Nodes are searched in M.nodes and <lNodesToValidNodes> is used to convert to valid node indices
      % that way, if not found, <lNodesToTCmdNodes> can be used to get the node index in the command vector
      [~, ind1] = ismember(M.Gdyn{i,1}, M.nodes);       % ind1=0 means the node is in the command vector
      % ismember supports char but cell array of char as well
      [~, ind2s] = ismember(M.Gdyn{i,2}, M.nodes);

      lValid_ind1 = lNodesToValidNodes(ind1);
      if lValid_ind1 == 0 % No more in valid node list
        lValid_ind1 = lNodesToTCmdNodes(ind1);
        assert(lValid_ind1 != 0);
        Gdyn(i,1:3) = {0, ind2s, lValid_ind1};
      elseif any(ind2s == 0)
        lValid_ind2s = lNodesToTCmdNodes(ind2s);
        assert(all(lValid_ind2s != 0));
        Gdyn(i,1:3) = {ind1, ind2s, lValid_ind2s};
      else
        Gdyn(i,1:3) = {ind1, ind2s, 0};
      endif

      Gdyn(i,4:5) = { M.Gdyn{i,3}, M.Gdyn{i,4} };
    endfor
  clear ind1 ind2 lValid_ind1 lValid_ind2s lNodesToValidNodes lNodesToTCmdNodes;
  endif

  % Resolution numérique

  userData = [];
  iT = 1;
  lTimeTrace = struct(  'updateMatrix', 0, ...
                        'solveMatrix', 0, ...
                        'prepareMatrix', 0, ...
                        'iterTime', [Inf 0 0 0]);  % min, max, mean, moving average

  % Set the first temperature based on initial temperature
  Tmat(:,1) = T0 + lConvertUnitOffset;
  clear T0;

  OptOptions = optimset('AutoScaling', 'off', ...
                        'MaxIter', 150, ...
                        'MaxFunEvals', 800*nNodes, ...
                        'Display', 'iter', ...
                        'TolX', 1E-8);

  % Init command vector and init command matrix
  % The base command matrix B is retrieve by GetCmdVEctor since for some command type
  % B is not changed after (only for k=1)
  if is_function_handle(options.iterstartfunc)
    userData = options.iterstartfunc(1,tStart,Tmat(:,1),[]);
  endif
  [Unow Bnow] = INT_GetCmdVector(Cmd, Umat(:,1), B, 1, tStart, userData, options); % Cmd, index, time
  [Gnow Bnow] = INT_UpdateMatrix(G, Unow, Bnow, Gdyn, Tmat(:,1), tStart, 1); % Init Gnow and Bnow matrix

  Umat(:,1) = Unow;

  lMatChanged = struct( 'G', true, ...
                        'Glast', true,...
                        'C', true, ...
                        'B', true, ...
                        'Blast', true, ...
                        'A', true, ...
                        'dt', true);

  if iscell(tVec)
    dt = tVec{2};
    lMatChanged.dt = false; % Always false in that case, and it will remain like this
  else
    dt = tVec(2)-tVec(1);
  endif

  ABackup = [];
  t = zeros(nt,1);
  t(1) = tStart;

  i=1;
  lIterationCount = 0;
  lFirstIteration = true;
  lNewTimeStep = false;     % Holds whether or not the current loop is an additionnal time step (not specified in the time vector)
  tCurrent = tStart;

  while i < nt
    lIterationCount++;
    lIterCpuTime = cputime();

    % If <tVec> is a cell, it means <dt> is constant
    if ~iscell(tVec)
      dtold = dt;   %tVec(i) - tVec(i-1);
      dt = tVec(i+1) - tCurrent;
      if ~isempty(options.maxdt) && (dt > 1.1*options.maxdt)
        dt = options.maxdt;
        tCurrent += dt;
        lNewTimeStep = true;
      else
        i++;
        tCurrent = tVec(i);
        lNewTimeStep = false;
      endif

      lMatChanged.dt = ( abs(dt - dtold)/(tEnd-tStart) > 1E-14 ) || lFirstIteration;
    else
      i++;
      tCurrent += dt;
    endif

    t(i) = tCurrent;

    Glast = Gnow;     lMatChanged.Glast = lMatChanged.G;
    Blast = Bnow;     lMatChanged.Blast = lMatChanged.B;
    Tlast = Tmat(:,iT);
    Ulast = Unow;

    % Call user function if specified
    if is_function_handle(options.iterstartfunc), userData = options.iterstartfunc(i,tCurrent,Tlast); endif;

    % Get command vector and reinit command matrix from B (flux type command)
    tmp = cputime();
    [Unow, Bnow, lMatChanged.B] = INT_GetCmdVector(Cmd, Unow, B, i, tCurrent, userData, options);
    [Gnow, Bnow, lMatChanged.G, lMatChanged.B] = INT_UpdateMatrix(G, Unow, Bnow, Gdyn, Tlast, tCurrent, i); % (temperature type command)
    lMatChanged.C = false;

    lTimeTrace.updateMatrix += cputime()-tmp;

    if lMatChanged.G || lMatChanged.dt || lFirstIteration || lMatChanged.C
      A = (C - Gnow*(dt/2));
      lMatChanged.A = true;
    else
      lMatChanged.A = false;
    endif

    if strcmpi(options.algorithm, 'exact')
      % Solve Ax=b
      b = (C + Glast*(dt/2))*Tlast + (dt/2)*(Blast + Bnow);

      % Solve min (residuals)^2
      if ~isempty(M.rad)
        fobj = @(T) sumsq(A*T - ...
                          b - ...
                          (dt/2)*RradTotal * ((T/273.15).^4 + (Tlast/273.15).^4) - ...
                          (dt/2)*RradU * ((Unow(BTemperatureIndex)/273.15).^4 + (Ulast(BTemperatureIndex)/273.15).^4));
      else
        fobj = @(T) sumsq(A*T - b);
      endif

      tmp = cputime();
      [X fval info] = fminunc(fobj, Tlast, OptOptions);
      lTimeTrace.solveMatrix += cputime()-tmp;
      disp(info);

    else

      % Apply radiation problem
      if ~isempty(M.rad)
        tmp = cputime();
        [dARad, dbRad] = INT_ApplyRadiationMatrix(lRadiationMatrix, lRadiationMatrixCmd, Tlast, Ulast, Unow, lNodeTempRadiationList);
        dARad *= dt;
        dbRad *= dt;

        ABackup = A;
        A += dARad;
        dARad = [];
        lMatChanged.A = true;
        lTimeTrace.prepareMatrix += cputime()-tmp;
      endif

      if strcmpi(options.algorithm, 'linear') % Solve A.T=b
        % Solve Ax=b
        tmp = cputime();
        b = (C + Glast*(dt/2))*Tlast + (dt/2)*(Blast + Bnow);

        if ~isempty(M.rad)
          b += dbRad;
        endif

  ##        ABackup = A;
  ##        lTLastExp1 = Tlast/273.15;
  ##        lTLastExp3 = lTLastExp1.^3;
  ##        A -= 2*dt*lRadiationMatrix * spdiags(lTLastExp3/273.15, 0, numel(lTLastExp3), numel(lTLastExp3)); % .* (Tlast'/273.15).^3/273.15; % Columns are multiplied by Tlast
  ##        b += -dt*lRadiationMatrix * (lTLastExp3 .* lTLastExp1);
  ##
  ##        if ~isempty(lRadiationMatrixCmd)
  ##          lTLastExp1= Ulast(lNodeTempRadiationList)/273.15;
  ##          lTLastExp3 = lTLastExp1.^3;
  ##          lTLastExp4 = lTLastExp3.*lTLastExp1;
  ##
  ##          b += dt*( -(lRadiationMatrixCmd * lTLastExp4) ...
  ##                    +2*(lRadiationMatrixCmd * (lTLastExp3 .* (Unow(lNodeTempRadiationList)/273.15))));
  ##        endif
##        endif
        lTimeTrace.prepareMatrix += cputime()-tmp;

        tmp = cputime();
          X = linsolve(A, b, struct('POSDEF', true));
        lTimeTrace.solveMatrix += cputime()-tmp;

      elseif strcmpi(options.algorithm, 'lu')
        assert(isempty(M.rad), 'Not supported yet');

        tmp = cputime();
        if lMatChanged.A
          % Solve Ax=b
          [LU_L, LU_U, LU_P, LU_Q] = lu(A, 'vector');
          LU_Q(LU_Q) = 1:numel(LU_Q); % = sparse(LU_Q, 1:nNodes, 1, nNodes, nNodes);
        endif

        if lMatChanged.C || lMatChanged.Glast || lMatChanged.dt
          LU_Z = C + Glast*(dt/2);
        endif

        %b = (C + Glast*dt/2)*Tlast + dt/2*(Blast*Ulast + Bnow*Unow);
        b = LU_Z*Tlast + (dt/2)*( Blast + Bnow );

        if ~isempty(M.rad), b += dbRad; endif

  ##      Matrix case
  ##      X = LU_L \ (LU_P*b);
  ##      X = LU_U \ X;
  ##      X = LU_Q * X;
  %     Vector case
        X = LU_L \ (b(LU_P));
        X = LU_U \ X;
        X = X(LU_Q);
        lTimeTrace.solveMatrix += cputime()-tmp;

      elseif strcmpi(options.algorithm, 'chol')
        assert(isempty(M.rad), 'Not supported yet');

        tmp = cputime();
        if lMatChanged.A
          tmp = cputime();

          % Solve Ax=b
          A = (C - Gnow*(dt/2));
          [CHOL_R, p, CHOL_Q] = chol(A, 'lower', 'vector');
          assert(p==0);

          % Revert the permutation vector, ie transpose Q.
          CHOL_QP = CHOL_Q;
          CHOL_QP(CHOL_Q) = 1:numel(CHOL_Q);

          CHOL_R = matrix_type(CHOL_R, "lower");
          CHOL_RT = CHOL_R';
          CHOL_RT = matrix_type(CHOL_RT, "upper");

          lTimeTrace.prepareMatrix += cputime()-tmp;
        endif

        if lMatChanged.C || lMatChanged.Glast || lMatChanged.dt
          CHOL_Z = C + Glast*(dt/2);
        endif

        b = CHOL_Z*Tlast + (dt/2)*(Blast + Bnow);

        if ~isempty(M.rad), b += dbRad; endif

        % Following commented lines are referring to the matrix version
        % when input parameter 'vector' is not used in chol function.
        % At present time, there is an error in Octave documentation.
        % In chol documentation, it is written: R' * R = Q' * A * Q.
        % but it is in fact: R * R' = Q' * A * Q.
  ##      X = CHOL_R \ (CHOL_Q'*b);
  ##      X = CHOL_R' \ X;
  ##      X = CHOL_Q * X;
        % With argument 'vector', previous sparsem matrix multiplication becomes
        % vector permutation
        X = CHOL_R \ (b(CHOL_Q));
        X = CHOL_RT \ X;
        X = X(CHOL_QP);
        lTimeTrace.solveMatrix += cputime()-tmp;
      endif

      if ~isempty(M.rad), A = ABackup; ABackup = []; endif

    endif


    tmp = cputime() - lIterCpuTime;
    lTimeTrace.iterTime(1) = min(lTimeTrace.iterTime(1), tmp);
    lTimeTrace.iterTime(2) = max(lTimeTrace.iterTime(2), tmp);
    lTimeTrace.iterTime([3 4]) += tmp;

    if strcmpi(options.info, 'iter')
      disp(sprintf('Iteration %d/%d done: %.1f ms ; total time %.1 s', i,nt, tmp*1000, toc(lRunningTimeId)));
    elseif (strcmpi(options.info, 'default') && (mod(i, 100) == 0))
      disp(sprintf('Iteration %d/%d done: mean iter time = %.1f ms ; total time: %.1f s', ...
              i,nt, lTimeTrace.iterTime(4)/100*1000, toc(lRunningTimeId)));
      lTimeTrace.iterTime(4) = 0;
    endif

    % Before modifying matrix T, the reference counter of matrix T must be decremented
    % to avoid that "T(:,iT) = X" forces a copy of the entire matrix since Tlast
    % is still referencing matrix T.
    Tlast = [];

    % Sometimes, if the time step computed from the time vector "t" specified by the user
    % is too large, it is cut to match <options.maxdt>, it that case, the current loop
    % must not be saved, otherwise the original matrix T will be too small.
    % Instead of resizing the matrix, it is chosen to skip the time step.
    if ~lNewTimeStep
      if options.all, iT++; endif;
      Tmat(:,iT) = X;

      if options.mergeTU, Umat(:,iT) = Unow; endif;
    endif

    lFirstIteration = false;

  endwhile


  if options.mergeTU
    Tmat = [Tmat; Umat];
    Gnodes = [Gnodes; Unodes];
  endif

  if lConvertToKelvin, Tmat -= lConvertUnitOffset; endif;

  M.timer.solvingTime += cputime() - lRunningTime;

  if options.verbose
    disp(sprintf('Solving model <%s> done in %.2f ms (%d iterations)', M.name, 1000*M.timer.solvingTime, lIterationCount));
    disp(sprintf('Stats: mean iteration time = %.1f ms, min = %.1f ms, max = %.1f ms', 1000*lTimeTrace.iterTime(3)/lIterationCount, 1000*lTimeTrace.iterTime(1), 1000*lTimeTrace.iterTime(2)));
    disp(sprintf('Update matrix time = %.2f ms', 1000*lTimeTrace.updateMatrix));
    disp(sprintf('Solve matrix time = %.2f ms', 1000*lTimeTrace.solveMatrix));
    disp(sprintf('Prepare matrix time = %.2f ms', 1000*lTimeTrace.prepareMatrix));
  endif;

endfunction










function RadiationMatrix = INT_PrepareRadiationMatrix(radModelList, nodes)
  assert(all(HT_CheckType(radModelList, 'radiationModel')));

  RadiationMatrix = sparse(numel(nodes), numel(nodes));

  for i=1:numel(radModelList)
    Mrad = radModelList{i};
    [~, ind] = ismember(Mrad.nodes, nodes);

    assert(iscolumn(Mrad.area) && iscolumn(Mrad.emissivity));

    Id = speye(size(Mrad.VF));
    R = (Mrad.area .* full(Mrad.VF - Id)) * ...
            inv(Id - spdiags(1 - Mrad.emissivity, 0, rows(Mrad.VF), columns(Mrad.VF)) * Mrad.VF) * ...
            (diag(Mrad.emissivity) * (5.67E-8 * 273.15^4));

    RadiationMatrix(ind,ind) += R;
  endfor
endfunction

function [dA, db] = INT_ApplyRadiationMatrix(lRadiationMatrix, lRadiationMatrixCmd, Tlast, Ulast, Unow, lNodeTempRadiationList)

  lTLastExp1 = Tlast/273.15;
  lTLastExp3 = lTLastExp1.^3;
##  dA = -2*lRadiationMatrix * spdiags(lTLastExp3/273.15, 0, numel(lTLastExp3), numel(lTLastExp3)); % .* (Tlast'/273.15).^3/273.15; % Columns are multiplied by Tlast
##  db = -lRadiationMatrix * (lTLastExp3 .* lTLastExp1);
##
##  if ~isempty(lRadiationMatrixCmd)
##    lTLastExp1= Ulast(lNodeTempRadiationList)/273.15;
##    lTLastExp3 = lTLastExp1.^3;
##    lTLastExp4 = lTLastExp3.*lTLastExp1;
##
##    db += ( -(lRadiationMatrixCmd * lTLastExp4) ...
##              +2*(lRadiationMatrixCmd * (lTLastExp3 .* (Unow(lNodeTempRadiationList)/273.15))));
##  endif
%09/06/2023 it seems there is an error
##  dA = -2*lRadiationMatrix * spdiags(lTLastExp3/273.15, 0, numel(lTLastExp3), numel(lTLastExp3)); % .* (Tlast'/273.15).^3/273.15; % Columns are multiplied by Tlast
##  db = -lRadiationMatrix * (lTLastExp3 .* lTLastExp1);
% Je ne comprends pas pour la linéarisation pose problème ici.
  dA = 0;
  db = lRadiationMatrix * (lTLastExp3 .* lTLastExp1);

  if ~isempty(lRadiationMatrixCmd)
    lULastExp1= Ulast(lNodeTempRadiationList)/273.15;
##    lULastExp3 = lULastExp1.^3;
##    lULastExp4 = lULastExp3.*lULastExp1;

    db += lRadiationMatrixCmd * (lULastExp1.^4 + (Unow(lNodeTempRadiationList)/273.15).^4)/2;
##    db += lRadiationMatrixCmd * (2*lULastExp1.^4 + 4*lULastExp1.^3.*(Unow(lNodeTempRadiationList)/273.15 - lULastExp1))/2;
  endif

endfunction



function [G B Gchanged Bchanged] = INT_UpdateMatrix(G, Unow, B, Gdyn, T, t, tindex)
  nDyn = size(Gdyn, 1);
  Gchanged = (nDyn > 0);
  Bchanged = false;

  for i=1:nDyn
    assert(any(Gdyn{i,5} == [0 2 3]), 'Invalid function type');
    m = Gdyn{i,1};
    ns = Gdyn{i,2};
    ks = Gdyn{i,3}; % ks is the index in command vector

    nsValid = (ns != 0);
    nsIndex = ns(nsValid);
    nsInvalid = find(~nsValid);

    % Evaluation of conductance
    if Gdyn{i,5} == 0 % Matrix case
      lConductanceVec = Gdyn{i,4}(:, max(tindex, size(Gdyn{i,4},2)));
    elseif Gdyn{i,5} == 2
      lConductanceVec = Gdyn{i,4}(t, tindex); % 2 two argument (time, index)
    else
      error('Not implemented');
      % If the function has 3 input arguments, T, t, i
      % One must first build the Temperature vector T
      % The first elt of T is the first node temperature (node "m")
      % The other elts are the node temperatures (node list "ns" where empty are taken into "ks")
      % Indeed, some node indices of "ns" could be empty, meaning that the node temperature is in fact command temperature
##      val = Gdyn{i,4}(Tparam, t, tindex);
##      assert(numel(val) == numel(ns), 'Conduction function must return a vector of size equal to the number of nodes');
    endif

    % If ns has several element, val is a vector of size ns, otherwise it is a scalar
    lValidConductance = lConductanceVec(nsValid);

    lDiag = spdiags(G, 0);
    lDiag(nsIndex) -= lValidConductance;
    G = spdiags(lDiag, 0, G);   % Apply the changes to the conductance matrix

    if m != 0 % If the first node is valid
      G(m,m)        -= sum(lConductanceVec);    % Sum over valid nodes and nodes inside command vector U
      G(nsIndex, m) += lValidConductance;
      G(m,nsIndex)  += lValidConductance;
      B(m)          += sum(lConductanceVec(nsInvalid) .* Unow(ks(nsInvalid)));  % Connect to node in the command vector U
      Bchanged = Bchanged || ~isempty(nsInvalid);
    else
      % m==0 means all nodes of ns are valid
      % Here, ks contains the index of the temperature in Unow vector
      B( nsIndex) += lValidConductance * Unow(ks);
      Bchanged = true;
    endif
  endfor
endfunction

% Returns the command vector U at time t,k and adjust command matrix (for flux type cmd)
% B: initial command matrix (not last command matrix)
% k: iteration index
% t: current time
function [Uvec B Bchanged] = INT_GetCmdVector(Cmd, Uvec, B, k, t, userData, options)
  HT_ImportConstants();

  % First iteration ?
  if k <= 1
    Bchanged = true;
  else
    Bchanged = false;
  endif

  % Command vector BT is similar to vector B, but it the temperature of each nodes
  % and checks that no command attempt to set severals times the same node. If it occurs, it checks that the temperature is
  % equal.
  Uref = NA(size(Uvec)); % Keep all node temperatures
  UrefSet = false(size(Uvec));

  for i=1:numel(Cmd)
    cmd = Cmd{i};

    Uloc = NA(size(Uvec));    % Cleared for each command and compared to BTref after a command is applied to B

    if cmd.sys.faceProvided
      for m=1:numel(cmd.nodes)

        switch cmd.sys.type(m)
          case HT_BTYPE_T_CONSTANT_UNIFORM        % Set during solveModel initialisation
          case HT_BTYPE_T_CONSTANT_NON_UNIFORM    % Set during solveModel initialisation
          case HT_BTYPE_T_VARIABLE_UNIFORM        %
            B += cmd.sys.BT{m} * cmd.data{m}(k);
            Uloc(cmd.sys.faceRemoved{m}) = cmd.data{m}(k);
            Bchanged = true;
            if options.mergeTU, Uvec( cmd.sys.faceRemoved{m} ) = cmd.data{m}(k); endif;

          case HT_BTYPE_T_VARIABLE_NON_UNIFORM    %
            lTmp = accumdim(cmd.sys.faceIndexOperation{m}, [cmd.data{m}(:,k) .* cmd.area{m} , cmd.area{m}]);
            lTmp = lTmp(:,1) ./ lTmp(:,2);
            B += cmd.sys.BT{k} * lTmp;
            Uloc(cmd.sys.faceRemoved{m}) = lTmp;
            Bchanged = true;
            if options.mergeTU, Uvec( cmd.sys.faceRemoved{m} ) = lTmp; endif;

          case HT_BTYPE_T_INTERPOLATION_ARRAY
            lTmp = interp1(cmd.data{m}(:,2), cmd.data{m}(:,1), t);
            B += cmd.sys.BT{m} * lTmp;
            Uloc(cmd.sys.faceRemoved{m}) = lTmp;
            Bchanged = true;
            if options.mergeTU, Uvec( cmd.sys.faceRemoved{m} ) = lTmp; endif;

          % ======= Flux type =========
          case HT_BTYPE_FLUX_CONSTANT_UNIFORM         % .sys.type=21 : flux constant and uniform (Set during solveModel initialisation)
          case HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM     % .sys.type=22 : variable and uniform flux (Set during solveModel initialisation)
          case HT_BTYPE_FLUX_VARIABLE_UNIFORM         % .sys.type=23 : constant and non uniform flux
            if isempty(cmd.area) || isempty(cmd.area{m})
##              B( cmd.sys.faceIndex{m} ) += cmd.data{m}(k);
              B += accumarray(cmd.sys.faceIndex{m}, cmd.data{m}(k), size(B));
            else
##              B( cmd.sys.faceIndex{m} ) += cmd.data{m}(k) .* cmd.area{m};
              B += accumarray(cmd.sys.faceIndex{m}, cmd.data{m}(k) .* cmd.area{m}, size(B));
            endif
            Bchanged = true;

          case HT_BTYPE_FLUX_VARIABLE_NON_UNIFORM     % .sys.type=24 : variable and non uniform flux
            if isempty(cmd.area) || isempty(cmd.area{m})
              B += accumarray(cmd.sys.faceIndex{m}, cmd.data{m}(:,k), size(B));
            else
              B += accumarray(cmd.sys.faceIndex{m}, cmd.data{m}(:,k) .* cmd.area{m}, size(B));
            endif
            Bchanged = true;
          otherwise
            error(sprintf('Unknown command type %d', cmd.sys.type));
        endswitch

        % This additionnal check is not done for flux command type since there
        % is no issue when several occurences of the same node exist
        if cmd.sys.type < HT_BTYPE_FLUX_CONSTANT_UNIFORM  % Temperature type command ?   ~isempty(cmd.sys.faceIndexUnique)
          lTemp = UrefSet(cmd.sys.faceRemoved{m});
          assert(all(~lTemp) || all(abs(Uloc(lTemp)-Uref(lTemp)) < eps('double')), 'Temperature command: Conflict of node temperature');
          Uref(lTemp) = Uloc(lTemp);
          UrefSet(lTemp) = true;
        endif
      endfor

    else % if cmd.sys.faceProvided

      switch cmd.sys.type
        case HT_BTYPE_T_CONSTANT_UNIFORM        % Set during solveModel initialisation
        case HT_BTYPE_T_CONSTANT_NON_UNIFORM    % Set during solveModel initialisation
        case HT_BTYPE_T_VARIABLE_UNIFORM        %
          B += cmd.sys.BT * cmd.sys.data(k);
          Uloc(cmd.sys.nodeRemoved) = cmd.sys.data(k);
          Bchanged = true;
          if options.mergeTU, Uvec( cmd.sys.nodeRemoved ) = cmd.data(k); endif;

        case HT_BTYPE_T_VARIABLE_NON_UNIFORM    %
          B += cmd.sys.BT * cmd.data(:,k);
          Uloc(cmd.sys.nodeRemoved) = cmd.sys.data(:,k);
          Bchanged = true;
          if options.mergeTU, Uvec( cmd.sys.nodeRemoved ) = cmd.data(:,k); endif;

        case HT_BTYPE_T_INTERPOLATION_ARRAY
          lTmp = interp1(cmd.data(:,2), cmd.data(:,1), t);
          B += cmd.sys.BT{m} * lTmp;
          Uloc(cmd.sys.nodeRemoved) = lTmp;
          Bchanged = true;
          if options.mergeTU, Uvec( cmd.sys.nodeRemoved ) = lTmp; endif;

        % ======= Flux type =========
        case HT_BTYPE_FLUX_CONSTANT_UNIFORM         % .sys.type=21 : flux constant and uniform (Set during solveModel initialisation)
        case HT_BTYPE_FLUX_CONSTANT_NON_UNIFORM     % .sys.type=22 : variable and uniform flux (Set during solveModel initialisation)
        case HT_BTYPE_FLUX_VARIABLE_UNIFORM         % .sys.type=23 : constant and non uniform flux
          if isempty(cmd.area)
            B(cmd.sys.nodeIndex) += cmd.data(k);
          else
            B(cmd.sys.nodeIndex) += cmd.data(k) .* cmd.area;
          endif
          Bchanged = true;

        case HT_BTYPE_FLUX_VARIABLE_NON_UNIFORM     % .sys.type=24 : variable and non uniform flux
          if isempty(cmd.area)
            B(cmd.sys.nodeIndex) += cmd.data(:,k);
          else
            B(cmd.sys.nodeIndex) += cmd.data(:,k) .* cmd.area;
          endif
          Bchanged = true;

        otherwise
          error(sprintf('Unknown command type %d', cmd.sys.type));
      endswitch

      % This additionnal check is not done for flux command type since there
      % is no issue when several occurences of the same node exist
      if cmd.sys.type < HT_BTYPE_FLUX_CONSTANT_UNIFORM  % Temperature type command ?   ~isempty(cmd.sys.faceIndexUnique)
        lTemp = UrefSet(cmd.sys.nodeRemoved);
        assert(all(~lTemp) || all(abs(Uloc(lTemp)-Uref(lTemp)) < eps('double')), 'Temperature command: Conflict of node temperature');
        Uref(lTemp) = Uloc(lTemp);
        UrefSet(lTemp) = true;
      endif

    endif % if cmd.sys.faceProvided

  endfor

endfunction


% Returns an array by merging A and B
% Y(i) = A(i) if not NAN, else B(i)
% Error if A(i) and B(i) not NaN and different
% indAtoB is an array of index that link nodes of A to nodes of B
% indAtoB(i) is the index in B of the node i (0 if it does not exist in B)
function A = MergeArray(A, B, indAtoB)
  assert(numel(A) == numel(indAtoB));

  AIsNan = isnan(A);
  ANotIsNan = ~AIsNan;          % Logical array
  ANotIsNan = find(ANotIsNan);  % Index array
  ANotIsNan(indAtoB(ANotIsNan) == 0) = []; % Remove nodes that does not exist in B

  AValue = A(ANotIsNan);
  BValue = B(indAtoB(ANotIsNan));
  assert(all( (AValue == BValue) | isnan(AValue) | isnan(BValue)), 'Some values of initial temperature do not match');

  AIsNan = find(AIsNan);        % Index array
  AIsNan(indAtoB(AIsNan) == 0) = [];
  A(AIsNan) = B(indAtoB(AIsNan));
endfunction

% Find the node name "node" in command nodes specified in Cmd list "Cmd"
% Only commands Cmd(Select) are taken into account.
function ind = Int_FindNodeInCmd(node, Cmd, Select)
  assert(nargin == 3);
  assert(ischar(node) || iscellstr(node));
  assert(HT_CheckType(Cmd, 'command'));

  if ischar(node), node = { node }; endif;

  ind = zeros(numel(node), 1); % Default return value
  if iscolumn(Select), Select = Select'; endif;

  for i=Select
    lInvalidIndex = find(ind == 0);
    if isempty(lInvalidIndex), break; endif;

    [tf, s_idx] = ismember(node(lInvalidIndex), Cmd{i}.nodes);
    tf(tf) = i;
    ind(lInvalidIndex) = tf;
  endfor
endfunction


