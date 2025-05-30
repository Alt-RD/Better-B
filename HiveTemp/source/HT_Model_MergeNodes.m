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
% Fusionne plusieurs noeuds entre eux et modifie les matrices
% Une fois la fusion réalisée, les noeuds fusionnés disparaissent et n'existeront
% plus dans le modèle. Ils seront remplacés par un unique noeud. Le nom de ce nouveau
% noeud est par défaut le nom du 1er noeud de l'ensemble fusionné. Il est possible
% de spécifier explicitement un nouveau nom.
%
% Arguments:
% .model = [model] l'objet du modèle dont on souhaite fusionner les noeuds.
% .nodeList = [cell of string] or [index vector] or [cell of (cell of string or index vector)]
%             Liste des noeuds à fusionner, spécifiés par nom ou par indice (plus rapide)
%             Si une liste contenant des listes de noeuds est spécifiée, alors plusieurs opérations
%             indépendantes de fusion sont réalisées. <newName> doit avoir la taille appropriée.
% .newName = [string] or [cell of string]
%
% Retourne une structure model modifié

function M = HT_Model_MergeNodes(M, nodeList, newName)
  assert(HT_CheckType(M, 'model'));

  % Test if nodeList is a cell array of string ?
  if iscell(nodeList)
    if all(cellfun(@(v) ischar(v), nodeList))
      nodeList = { nodeList };
      assert(numel(newName) == 1, 'Invalid size of <newName>');
    endif
  else
    assert(isinteger(nodeList) && (iscolumn(nodeList) || isrow(nodeList)), 'Invalid index vector');
    nodeList = { nodeList };
    assert(all(cellfun(@(v) ischar(v) && ~isempty(v), nodeList)), 'Invalid parameter <newName>. Valid node names must be specified');
  endif

  if ischar(newName), newName = { newName }; endif
  if isempty(newName)
    newName = cell(numel(nodeList), 1);
  endif

  assert(iscell(newName) && all(cellfun(@(v) ischar(v), newName)), 'Invalid type of <newName>. Must be cell array');
  assert(numel(newName) == numel(nodeList), 'Invalid size of <newName> and <nodeList>');

  nNodes = numel(M.nodes);
  lToBeRemoved = false(nNodes, 1);

  for i=1:numel(nodeList)
    if iscell(nodeList{i})
      [lFound, lIndexInModel] = ismember(nodeList{i}, M.nodes);
      assert(all(lFound), sprintf('Some nodes (%d) do not exist in the model <%s>', sum(~lFound), M.name));
      nodeList{i} = lIndexInModel;
    else
      assert(all(nodeList{i} > 0) && (max(nodeList{i}) <= nNodes), 'Some nodes index are invalid');
    endif
  endfor

  for i=1:numel(nodeList)
    lNodeSet = false(nNodes, 1);
    lNodeSet(nodeList{i}) = true;

    assert(sum(lNodeSet) == numel(nodeList{i}), 'Duplicate nodes were specified');
    lRefIndex = nodeList{i}(1);
    lNodeSet(lRefIndex) = false;

    lMergeLine =    sum(M.G(lNodeSet,:), 1);
    M.G(lRefIndex,:) += lMergeLine;
    M.G(:,lRefIndex) += lMergeLine';
    M.G(lRefIndex,lRefIndex) = 0;
    M.G(lRefIndex,lRefIndex) = -sum(M.G(lRefIndex,~lNodeSet));

    M.C(lRefIndex) += sum(M.C(lNodeSet));
    M.T0(lRefIndex) = mean([M.T0(lRefIndex); M.T0(lNodeSet)]);

    if ~isempty(newName{i})
      M.nodes{lRefIndex} = newName{i};
    endif

    % It seems necessary...
    M.G(lNodeSet,:) = 0;
    M.G(:,lNodeSet) = 0;

    lToBeRemoved |= lNodeSet;
  endfor

  M.G(lToBeRemoved,:) = [];
  M.G(:,lToBeRemoved) = [];

  M.T0(lToBeRemoved) = [];
  M.C(lToBeRemoved) = [];

  M.nodes(lToBeRemoved) = [];

endfunction
