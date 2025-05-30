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
% Merge two model structure <M> and <Madd>. This functions adds the submodel
% <Madd> to a main model <M>
%
% Input arguments:
% .M     = [struct] data structure of type <model>
% .Madd  = [struct] data structure of type <model>
function M = HT_Model_Merge(M, Madd, options)
  lTicId = tic();

  assert(nargin >= 2, 'There must be 2 input arguments');

  if nargin < 3, options = struct(); endif;
  options = HT_CheckField(options, 'verbose', false,  @(v) islogical(v));

  if options.verbose, disp(sprintf('Merging models <%s> and <%s>', M.name, Madd.name)); endif;

  assert(isfield(M, 'G'), 'Field <G> is missing in model M');
  assert(isfield(Madd, 'G'), 'Field <G> is missing in model Mnew');

  assert(isfield(M, 'C'), 'Field <C> is missing in model M');
  assert(isfield(Madd, 'C'), 'Field <C> is missing in model Mnew');

  assert(isfield(M, 'Gdyn'), 'Field <Gdyn> is missing in model M');
  assert(isfield(Madd, 'Gdyn'), 'Field <Gdyn> is missing in model Mnew');

  assert(isfield(M, 'nodes'), 'Field <nodes> is missing in model M');
  assert(isfield(Madd, 'nodes'), 'Field <nodes> is missing in model Mnew');

  if ~isfield(M, 'submodel'), M = setfield(M, 'submodel', {}); endif;

  if ~isfield(options, 'verbose'), options.verbose = false;
    else options.verbose = cast(options.verbose, 'logical'); endif;

  % Check that the model was not already inserted
  if ~isempty(M.submodel)
  assert(~any(cellfun(@(m) strcmpi(Madd.name,m.name), M.submodel)), ...
    sprintf('The model <%s> was already inserted', Madd.name));
  endif

  nt = numel(M.nodes); %size(M.G,1); % Nombre de noeuds (base model)
  ntadd = numel(Madd.nodes); %size(Madd.G,1); % Nombre de noeuds (to be added)
##  nt = numel(M.nodes,1); % Nombre de noeuds (base model)
##  ntadd = numel(Madd.nodes,1); % Nombre de noeuds (to be added)

  Gnode1to2 = zeros(nt, 1);
  Gnode2to1 = zeros(ntadd, 1);

  [tf, ia] = ismember(M.nodes, Madd.nodes);
  Gnode1to2(tf) = ia(tf);
  Gnode2to1(Gnode1to2(tf)) = find(tf);

  lNewNodes = find(Gnode2to1 == 0);
  Gnode2to1(lNewNodes) = (numel(Gnode1to2)+1):(numel(Gnode1to2) + numel(lNewNodes));
  Gnode1to2 = [Gnode1to2; lNewNodes];

  M.nodes = [M.nodes; Madd.nodes(lNewNodes)];

##  Gnode1to2 = zeros(nt, 1);
##  Gnode2to1 = zeros(ntadd, 1);
##  for i=1:nt
##    ind = strcmpi(M.nodes{i}, Madd.nodes);
##    assert(any(numel(nonzeros(ind)) == [0 1])); % A node can mot appear twice
##
##    ind = find(ind);
##
##    if ~isempty(ind)
##      Gnode1to2(i) = ind;
##      Gnode2to1(ind) = i;
##    endif
##  endfor
##
##  for i=1:ntadd
##    if all(Gnode1to2 != i) %% It means Madd.notes(i) does not exist in M
##      M.nodes = [M.nodes; Madd.nodes(i)]; % On l'ajoute
##      Gnode1to2 = [Gnode1to2; i];
##      Gnode2to1(i) = numel(M.nodes);
##    endif
##  endfor
##
##  if (numel(Gnode1to2) != numel(Gnode1to2test)) || ~all((Gnode1to2 == Gnode1to2test) && (Gnode2to1 == Gnode2to1test))
##    disp('ici');
##  endif

% Unoptimized version
##  for i=1:ntadd
##    if all(Gnode1to2 != i) %% It means Madd.notes(i) does not exist in M
##      M.nodes = [M.nodes; Madd.nodes(i)]; % On l'ajoute
##      Gnode1to2(end+1) = i;
##      Gnode2to1(i) = numel(M.nodes);
##    endif
##  endfor

  ntnew = numel(M.nodes);  % New size with newly inserted nodes

  if options.verbose
    lCommonNodeCount = ntadd - numel(lNewNodes);
    if lCommonNodeCount > 0
      disp(sprintf('They are %d shared nodes between model <%s> and <%s>', lCommonNodeCount, M.name, Madd.name));
    endif;
  endif

  % Merge capacity matrix
  if isempty(M.C)
    M.C = NaN(ntnew, 1);
  else
    M.C = postpad(M.C, ntnew, NaN);
  endif

  % NaN values in M.C are replaced by values in Madd.C (maybe they are defined)
  M.C = MergeArray(M.C, Madd.C, Gnode1to2);

  % Merge conductance matrix
  if isempty(M.G)
    M.G = sparse(ntnew, ntnew);
  else
    M.G = resize(M.G, ntnew, ntnew);
  endif

  for i=1:size(Madd.G,1) % For each row of the conductance matrix
    newRow = Gnode2to1(i);
    M.G(newRow, Gnode2to1) += Madd.G(i,:);
  endfor

  % Merge time varying conductance
  M.Gdyn = [M.Gdyn; Madd.Gdyn];

  % Merge initial temperature vector
  if isempty(M.T0)
    M.T0 = NaN(ntnew, 1);
  else
    M.T0 = postpad(M.T0, ntnew, NaN);
  endif

  M.T0 = MergeArray(M.T0, Madd.T0, Gnode1to2);

  for i=1:numel(Madd.rad)
    assert(~any( cellfun(@(v) strcmpi(v.name, Madd.rad{i}.name), M.rad) ), 'Radiation problem is already added');
  endfor

  M.rad = [M.rad; Madd.rad]; % Add radiation problem
  M.submodel = [M.submodel; {Madd}];

  M.timer.buildingTime += Madd.timer.buildingTime;
  M.timer.solvingTime += Madd.timer.solvingTime;
  M.counter.warnings += Madd.counter.warnings;
  M.counter.errors += Madd.counter.errors;

  lExecTime = toc(lTicId);
  M.timer.buildingTime += lExecTime;

  if options.verbose, disp(sprintf('Merging models done in %.2f ms', 1000*lExecTime)); endif;

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
  assert(all( (AValue == BValue) | isnan(AValue) | isnan(BValue)), 'Some elements of the array can not be merged. Could be a conflict due to invalid node names.');

  AIsNan = find(AIsNan);        % Index array
  AIsNan(indAtoB(AIsNan) == 0) = [];
  A(AIsNan) = B(indAtoB(AIsNan));
endfunction

