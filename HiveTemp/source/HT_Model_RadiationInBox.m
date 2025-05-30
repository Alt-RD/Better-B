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

% Return a model of radiation heat transfer inside a box
% The box can not be oriented. Xm is the yz face in global axis
% Input arguments
% name: name of the model
% params: structure of parameters:
% .face = cell of face structures (order: Xm, Xp, Ym, Yp, Zm, Zp) in global axis
% .emissivity = vector dim 1x6 of emissivities of each face
%
% options: structure of options
% Output arguments:
% - M: Model object
% - F: Flux matrix
function [M F] = HT_Model_RadiationInBox(name, params, options)
  HT_ImportConstants();

  assert(nargin >= 2, 'There must be at least 2 input arguments');
  assert(numel(params) == 1, 'Invalid structure of parameters. numel > 1');

  if nargin < 3, options = struct(); endif;

  % Check input structure field names
  if HT_VAR_STRICT_INPUT_PARAMETERS
    lFieldNames = {'face', 'emissivity', 'wavelength'};
    lFieldTest = fieldnames(params);
    tf = ismember(lFieldTest, lFieldNames);
    assert(all(tf), sprintf('Invalid parameter name: %s', strjoin(cellfun(@(v) sprintf('%s',v), lFieldTest(tf), 'UniformOutput', false))));
    clear lFieldNames lFieldTest tf;
  endif

  % Set default structure parameter
  if HT_VAR_CHECK_INPUT
    options = HT_CheckField(options, 'verbose',       true  , @(v) islogical(v));
    options = HT_CheckField(options, 'verboseLevel',  0     , @(v) isscalar(v));

    params = HT_CheckField(params, 'face',        [], {"exist", @(v) (numel(v) == 6)});
    params = HT_CheckField(params, 'emissivity',  [], {@(v) (numel(v) == 6) && all(v > 0) && all(v <= 1)});
    params = HT_CheckField(params, 'wavelength',  [], {@(v) isscalar(v) && (v > 0.1) && (v < 100.0)}); % Unit is um

    if isstruct(params.face)
      params.face = arrayfun(@(v) v, params.face, 'UniformOutput', false);
    endif

    params.face = params.face(:);

    assert(abs(params.face{1}.norm + params.face{2}.norm) < 1E-10, 'Invalid faces XM/XP');
    assert(abs(params.face{3}.norm + params.face{4}.norm) < 1E-10, 'Invalid faces YM/YP');
    assert(abs(params.face{5}.norm + params.face{6}.norm) < 1E-10, 'Invalid faces ZM/ZP');
  endif

  if ~isempty(params.wavelength)
    params.emissivity = cellfun(@(v) HT_Material_GetEmissivity(v.material, params.wavelength), params.face);
  endif

  assert((numel(params.emissivity) == 6) && all(params.emissivity > 0) && all(params.emissivity <= 1), 'Invalid emissivity');

  % Get the total node count
  nNodesVec = cellfun(@(v) numel(v.nodes), params.face );
  nNodesVec = nNodesVec(:);
  assert(all(nNodesVec > 0), 'Some faces have no mesh');

  nNodesCum = cumsum([0; nNodesVec]); % Cumulative sum for index range computation
  nodesRange = [nNodesCum(1:(end-1))+1, nNodesCum(2:end)];
  nNodes = nNodesCum(end); % Total node count

  M = HT_Model_Init("radiationBox", name, params);
  Mrad = HT_Model_Radiation_Init(name);
  Mrad.VF = zeros(nNodes, nNodes); % Matrix of shape view factor

  % Create vectors containing node names, node area, node emissivity
  Mrad.nodes = cell2mat(cellfun(@(v) v.nodes, params.face, 'UniformOutput', false));
  Mrad.area = cell2mat(HT_Face_GetAbsoluteData(params.face, 'nodesArea'));
  Mrad.emissivity = repelem(params.emissivity, nNodesVec);

  % Adjust local axis of each face to make sure they match with the global axis
  % It allows distance between nodes to be computed.
  lAxis = eye(3);
  params.face{FaceXM} = HT_Face_AdjustAxis(params.face{FaceXM}, lAxis(:,[Y Z]));
  params.face{FaceXP} = HT_Face_AdjustAxis(params.face{FaceXP}, lAxis(:,[Y Z]));
  params.face{FaceYM} = HT_Face_AdjustAxis(params.face{FaceYM}, lAxis(:,[X Z]));
  params.face{FaceYP} = HT_Face_AdjustAxis(params.face{FaceYP}, lAxis(:,[X Z]));
  params.face{FaceZM} = HT_Face_AdjustAxis(params.face{FaceZM}, lAxis(:,[X Y]));
  params.face{FaceZP} = HT_Face_AdjustAxis(params.face{FaceZP}, lAxis(:,[X Y]));

  assert(all(params.face{FaceXM}.size == params.face{FaceXP}.size));
  assert(all(params.face{FaceYM}.size == params.face{FaceYP}.size));
  assert(all(params.face{FaceZM}.size == params.face{FaceZP}.size));

  lBoxDim = [params.face{FaceYM}.size(1) ; ...
             params.face{FaceXM}.size(1) ; ...
             params.face{FaceXM}.size(2)];
  lBoxArea = [prod(lBoxDim([2 3])) prod(lBoxDim([1 3])) prod(lBoxDim([1 2]))];

  assert(all(abs(params.face{FaceXM}.size - lBoxDim([2 3])) < 1E-14));
  assert(all(abs(params.face{FaceXP}.size - lBoxDim([2 3])) < 1E-14));
  assert(all(abs(params.face{FaceYM}.size - lBoxDim([1 3])) < 1E-14));
  assert(all(abs(params.face{FaceYP}.size - lBoxDim([1 3])) < 1E-14));
  assert(all(abs(params.face{FaceZM}.size - lBoxDim([1 2])) < 1E-14));
  assert(all(abs(params.face{FaceZP}.size - lBoxDim([1 2])) < 1E-14));

  lNodesDims = cell(6,1); % Avoid to perform the same multiplication severals times
  lNormVec = [ 1  -1  2  -2   3  -3];
  lComVec = zeros(6,6);  % xm xp ym yp zm zp
  lComVec(1,[3 4 5 6]) = [3 3 2 2];
  lComVec(2,[3 4 5 6]) = [3 3 2 2];
  lComVec(3,[5 6]) = [1 1];
  lComVec(4,[5 6]) = [1 1];

  for i=int32(1:6)
    lNodesDims{i} = params.face{i}.dims .* params.face{i}.size([1 1 2 2])';
    indDir = idivide(i-1,2, 'floor')+1; % Index of direction 1 for XM/XP, 2 for YM/YP, 3 for ZM/ZP

    for j=int32((i+1):6)
      if i == 1 % Initialize lNodesDims during the first iteration of i
        lNodesDims{j} = params.face{j}.dims .* params.face{j}.size([1 1 2 2])'; %repelem(params.face(j).size', 2);
      endif

      % Radiation view factor between face i and face j
      if any(i == [1 3 5]) && (j == i+1)
        % Faces are face to face
        % f1: Positions du rectangle 1
        %   > f1 = xm xp ym yp
        n1 = size(lNodesDims{i},1);
        n2 = size(lNodesDims{j},1);
        ldims1 = repelem(lNodesDims{i}, n2, 1);
        ldims2 = repmat(lNodesDims{j}, n1, 1);

        F12 = HT_ViewFactor_Rect2RectFace(ldims1, ldims2, lBoxDim(indDir));
        F12 = reshape(F12, n2, n1); % (F12)i,j = Fj,i (face 1 to face 2)
        % Chaque colonne de F12 correspond à un node différent de face 1 et
        %   balaie toutes les nodes de face 2
        % Compute surface area of each nodes and compare to the box face area
        S1 = (lNodesDims{i}(:,2)-lNodesDims{i}(:,1)) .* (lNodesDims{i}(:,4)-lNodesDims{i}(:,3));
        S2 = (lNodesDims{j}(:,2)-lNodesDims{j}(:,1)) .* (lNodesDims{j}(:,4)-lNodesDims{j}(:,3));

        % Compute view factor of face 2 nodes to face 1 nodes
        % eq to  F21 = diag(1./S2) * F21 * diag(S1)
        F21 = (F12 ./ S2 ).*(S1');

        % Check by comparison with 2 rectangle face to face
        F1B = sum(F12, 1)'; % Compute view factor between node1(face1) to all nodes of face2
        SA = sum(S1);
        SAchk = lBoxArea(indDir);
        assert(abs(SA-SAchk) < 1E-10, sprintf('SA area check failed: %e != %e', SA, SAchk));

        % Compute surface area of each nodes (face B) and compare to the box face area
        SBchk = lBoxArea(idivide(j,2, 'floor'));
        SB = sum(S2);
        % Fi,B was computed before (F1B). FB,A = 1/SB * sum ( FB,i / Si )
        FBA = 1/SB * sum(S1 .* F1B);
        % SA*FA,B = SB*FB,A then:
        FAB = 1/SA * sum(S1 .* F1B);
        FABchk = HT_ViewFactor_RectangularPlates( lBoxDim(mod(indDir,3)+1), ...
                                                  lBoxDim(mod(indDir+1,3)+1), ...
                                                  lBoxDim(indDir));
        assert(abs(FABchk - FAB) < 1E-10);
        clear n1 n2 ldims1 ldims2 F1B S1 SA SAchk SBchk SB FBA FAB;
        Mrad.VF( (nodesRange(i,1):nodesRange(i,2)),(nodesRange(j,1):nodesRange(j,2)) ) = F12';
        Mrad.VF( (nodesRange(j,1):nodesRange(j,2)),(nodesRange(i,1):nodesRange(i,2)) ) = F21;

      else
        % xm xp ym yp zm zp
        luVec = lAxis(:,lComVec(i,j));
        lv1Vec = lAxis(:,abs(lNormVec(j))) * sign(lNormVec(j));
        lv2Vec = lAxis(:,abs(lNormVec(i))) * sign(lNormVec(i));

        f1 = HT_Face_AdjustAxis(params.face{i}, [luVec lv1Vec]);
        f2 = HT_Face_AdjustAxis(params.face{j}, [luVec lv2Vec]);

        f1.dims .*= f1.size([1 1 2 2])'; %repelem(f1.size', 2);
        f2.dims .*= f2.size([1 1 2 2])'; %repelem(f2.size', 2);

        n1 = rows(f1.dims);
        n2 = rows(f2.dims);
        % [3 4 1 2] shift the order of components
        % Distance to the common axis first, then position along the common axis
        ldims1 = repelem(f1.dims(:, [3 4 1 2]), n2, 1);
        ldims2 = repmat(f2.dims(:, [3 4 1 2]), n1, 1);

        F12 = HT_ViewFactor_Rect2RectAngle(ldims1, ldims2);
        F12 = reshape(F12, n2, n1);

        % Compute surface area of each nodes and compare to the box face area
        S1 = (f1.dims(:,2)-f1.dims(:,1)) .* (f1.dims(:,4)-f1.dims(:,3));
        S2 = (f2.dims(:,2)-f2.dims(:,1)) .* (f2.dims(:,4)-f2.dims(:,3));
        % Compute view factor of face 2 nodes to face 1 nodes
        % eq to  F21 = diag(1./S2) * F21 * diag(S1)
        F21 = (F12 ./ S2 ).*(S1');

        % Check by computing the view factor of the whole face1 to face2
        F1B = sum(F12, 1)'; % Compute view factor between node1(face1) to all nodes of face2
        SA = sum(S1);
        SAchk = lBoxArea(indDir);
        if (abs(SA-SAchk) > 1E-10)
          assert(abs(SA-SAchk) < 1E-10);
        endif

        % Compute surface area of each nodes (face B) and compare to the box face area
        SBchk = lBoxArea(idivide(j,2, 'floor'));
        SB = sum(S2);
        % Fi,B was computed before (F1B). FB,A = 1/SB * sum ( Si * Fi,B )
        FBA = 1/SB * sum(S1 .* F1B);
        % SA*FA,B = SB*FB,A then:
        FAB = SB/SA * FBA;
        % The long index of lBoxDim just pick the right component in box dimension vector
        % with respect to the faces
        FABchk = HT_ViewFactor_RectangularPlatesAngle(...
                      lBoxDim(find(abs(cross(lAxis(:,lComVec(i,j)), lAxis(:,abs(lNormVec(j))) * sign(lNormVec(j)) )) > 0.5)), ...
                      lBoxDim(find(abs(cross(lAxis(:,lComVec(i,j)), lAxis(:,abs(lNormVec(i))) * sign(lNormVec(i)) )) > 0.5)), ...
                      lBoxDim(lComVec(i,j)));
        assert(abs(FABchk - FAB) < 1E-10);

        clear n1 n2 ldims1 ldims2 S1 S2 n1 n2 f1 f2 luVec lv1Vec lv2Vec FABchk FAB FBA SB SBchk F1B SA SAchk;
        Mrad.VF( (nodesRange(i,1):nodesRange(i,2)),(nodesRange(j,1):nodesRange(j,2)) ) = F12';
        Mrad.VF( (nodesRange(j,1):nodesRange(j,2)),(nodesRange(i,1):nodesRange(i,2)) ) = F21;
      endif
    endfor
  endfor

  M.rad = { Mrad }; % Add the radiation model to the model structure
  M.nodes = unique(Mrad.nodes, 'stable');
  M.G = []; % There are no conductance in this model
  M.C = NaN(numel(M.nodes), 1); % Capacity are not defined here
  M.T0 = NaN(numel(M.nodes), 1); % No initial temperature defined here

  % Check that the sum of every components on each row is equal to 1.
  assert(abs(sum(Mrad.VF,2)-1) < 1E-10, 'Some row have an invalid sum');

  % If the 2nd output argument is required
  % the flux matrix is computed and returned
  if isargout(2)
    Id = speye(size(Mrad.VF));
##    F = (Mrad.area .* full(Mrad.VF - Id)) * ...
##            (inv(Id - (1 - Mrad.emissivity) .* Mrad.VF) .* (Mrad.emissivity'));
    F = (full(Mrad.VF - Id)) * ...
            (inv(Id - (1 - Mrad.emissivity) .* Mrad.VF) .* (Mrad.emissivity'));
  endif

endfunction


