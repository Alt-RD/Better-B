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
%  Copyright (c) 2022 Montpellier-University
%  Copyright (c) 2023-2025 AltRD-Emmanuel Ruffio
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
% Solve the 1D water vapour transfer equations based on
% [Asseko Ella & al, Joseph Gril, 2024, "Diffusion properties of Gabonese tropical
% hardwoods and European softwoods measured with low-tech equipment and by an
% inverse method", Bois et Forêts des Tropiques – ISSN : L-0006-579X
% Volume 360 – 2e trimestre – juin 2024 – p. 65-76]
%
% This model is specified as "standalone" since it does to use the same logic
% as the other models of the HiveTemp library. It is mainly used for test,
% comparison, and evaluation purpose.
%
% Input arguments:
% .params = [struct] data structure containing following fields
%         .length = [double] length of the wall [m]
%         .type = ['isotherm'] configuration
%                 'isotherm' refers to the diffusion coefficient inside the sample.
%         .diffusion = [function_handle D = @(t,i,W)]
%                     or [2x1 vector A] with D=A(1)*(1+A(2)*W)
%         .convection = [1x2] vector containing the constant convection coefficient
%                             on each sides of the sample
%                     or [nt x 2] matrix of convection coefficients
%                             for both sides at each time step
%                     or [nt x 3] (t, h1, h2) matrix of convection coefficients
%                             for both sides with respect to time
%         .weq    =  [1x2] water content at equilibrium for each sides of the sample
%                 or [nt x 2] water content at equilibrium for each time step
%                 or [nt x 3] (t, weq1, weq2) water content at equilibrium with
%                             respect to time
%         .winit  = [1x1] or [nx1] initial water content, uniform or space varying
%                                  respectivitly (with <n> the number of nodes)
%         .n      = [integer] number of nodes
%         .u      = (optional) [dim nx1] in [0;1] specify the node borders using relative position
%         .gridType = ['ff', 'hf', 'fh', 'hh'] by default 'ff'.
%                     For numerical stability, 'ff' should be preferred. Otherwise
%                     large convection conductance compared to low water diffusion
%                     may generate oscillations.
%         .t      = [dim nt x 1] Time vector
% .options = [struct] data structure containing following options:
%         .verbose       = [logical] [true] display information during the calculations
%         .verboseLevel  = [integer] increasing this value increases the amount of
%                          information printed into the consol.
%
% Output arguments:
% .Wmat = [matrix nNode x nt] containing the water content at any times and any position
% .infos = [struct] structure containing various information
%       .nodePosition: normalized node position [0;1] inside the 1D geometry
%       .nodeVolume: normalized node volume (multiply by <length> to get the absolute volume)
function [Wmat infos] = HT_Standalone_Solve_Hygro1D(lParams, lOptions)
  HT_ImportConstants();

  lTicId = tic();
  infos = struct( 'nodePosition', [], ...
                  'nodeVolume', []);

  % Set default parameters for structure <options>
  lOptions = HT_CheckField(lOptions, 'verbose',        true,   { @(v) islogical(v) });
  lOptions = HT_CheckField(lOptions, 'verboseLevel',   double(lOptions.verbose),      { @(v) isscalar(v) && (v >= 0) && (round(v) == v) } );
  lOptions = HT_CheckField(lOptions, 'skipFieldCheck', false,    @(v) islogical(v));
  if lOptions.verbose, disp(sprintf('Solving standalone model hygro 1D...')); endif;

  % Check parameter structure fields. It is desirable to check that no silly errors
  % are made on the parameter name
  if ~lOptions.skipFieldCheck
    lValidParamList = { 'length', 'type', 'diffusion', 'convection', 'weq', 'winit', 'n', 'u', 'gridType', 't' };
    lParamList = fieldnames(lParams);
    lValidParam = cellfun(@(v) any(strcmp(v, lValidParamList)), lParamList);
    assert(all(lValidParam), sprintf('Invalid parameter name detected <%s>', strjoin(lParamList(~lValidParam), ',')));
    clear lParamList lValidParam lValidParamList;
  endif

  % Set default parameter values
  lParams = HT_CheckField(lParams, 'n',              [],                 @(v) isnumeric(v) && isscalar(v) && (round(v) == v) && v > 0);
  lParams = HT_CheckField(lParams, 'u',              [],                 @(v) isnumeric(u) && iscolumn(u) && all(u >= 0.0 & u <= 1.0) && all(diff(u) > 0) && (numel(u) == lParams.n || isempty(lParams.n)));
  if isempty(lParams.n), lParams.n = max(numel(lParams.u)-1, 0); endif;

  lParams = HT_CheckField(lParams, 'length',         [],                 {"exist", @(v) isnumeric(v) && isscalar(v) });
  lParams = HT_CheckField(lParams, 'type',           'isotherm',         {@(v) ischar(v) && any(strcmpi(v, 'isotherm'))});
  lParams = HT_CheckField(lParams, 'diffusion',      [],                 {"exist", @(v) is_function_handle(v) || (all(isnumeric(v(:))) && ((isrow(v) && numel(v)==2) || (columns(v)==3)))  });
  lParams = HT_CheckField(lParams, 'convection',     [],                 {"exist", @(v) is_function_handle(v) || (all(isnumeric(v(:))) && ((isrow(v) && numel(v)==2) || (columns(v)==3))) });
  lParams = HT_CheckField(lParams, 'weq',            [],                 {"exist", @(v) all(isnumeric(v(:))), @(v) (isrow(v) && numel(v)==2) || (columns(v)==3)});
  lParams = HT_CheckField(lParams, 'winit',          [],                 {"exist", @(v) all(isnumeric(v(:))) && any(numel(v) == [1 lParams.n]) });
  lParams = HT_CheckField(lParams, 'gridType',       'ff',               @(v) ischar(v) && any(strcmpi(v, {'ff', 'hf', 'fh', 'hh'})));
  lParams = HT_CheckField(lParams, 't',              [],                 {'exist', @(v) iscolumn(v) && all(diff(v) > 0) });

  % Structure containing all information and working variable for the calculations
  lLayer = struct('u', [], ...
                  'upos', [], ...
                  'usize', [], ...
                  'winit', lParams.winit, ...
                  'x_before', [], ...
                  'x_next', [], ...
                  'geometry', [],...
                  'command', [], ...
                  'nodes', {{}}, ...
                  'diffusion', lParams.diffusion, ...
                  'convection_left', [], ...
                  'convection_right', [], ...
                  'weq_left', [], ...
                  'weq_right', [] ...
                  );
  lLayer.u = lParams.u;

  % Compute the nodes positions and boundaries
  % <u> contains the nodes boundaries
  if isempty(lLayer.u) % If not specified by the user, regularly spaced nodes is used
    ldx = 1 / (lParams.n - 0.5*sum(lParams.gridType == 'h'));
    lLayer.usize = repmat(ldx, lParams.n, 1);
    if (lParams.gridType(1) == 'h'), lLayer.usize(1) *= 0.5; endif
    if (lParams.gridType(2) == 'h'), lLayer.usize(end) *= 0.5; endif
    lLayer.u = [0; cumsum(lLayer.usize)];
    clear ldx;
  else
    lLayer.n = numel(lParams.u) - 1;
    lLayer.usize = lParams.u(2:end) - lParams.u(1:(end-1));
  endif

  lLayer.upos = 0.5*(lLayer.u(1:(end-1)) + lLayer.u(2:end));
  if (lParams.gridType(1) == 'h'), lLayer.upos(1) = 0; endif
  if (lParams.gridType(2) == 'h'), lLayer.upos(end) = 1; endif

  if isscalar(lLayer.winit), lLayer.winit = repmat(lLayer.winit, lParams.n, 1); endif

  lLayer.x_before = (lLayer.upos - lLayer.u(1:(end-1)))*lParams.length;
  lLayer.x_next = (lLayer.u(2:end) - lLayer.upos)*lParams.length;

  if isnumeric(lLayer.diffusion) && isvector(lLayer.diffusion)
    lLayer.invDiffusion = @(t,i,W) 1./(lParams.diffusion(1)*(1+lParams.diffusion(2)*W));
  else
    lLayer.invDiffusion = @(t,i,W) 1./lParams.diffusion(t,i,W);
  endif

  if isrow(lParams.convection)
    a = lParams.convection;
    lLayer.convection_left = @(t, i, W) 1./a(1);
    lLayer.convection_right = @(t, i, W) 1./a(2);
  elseif columns(lParams.convection) == 3
    hMat = lParams.convection;
    lLayer.convection_left = @(t, i, W) 1./interp1(hMat(:,1), hMat(:,2), t);
    lLayer.convection_right = @(t, i, W) 1./interp1(hMat(:,1), hMat(:,3), t);
  elseif columns(lParams.convection) == 2
    hMat = lParams.convection;
    lLayer.convection_left = @(t, i, W) 1./hMat(i,1);
    lLayer.convection_right = @(t, i, W) 1./hMat(i,2);
  else
    error('Invalid parameter <lParams.convection>');
  endif

  if isrow(lParams.weq)
    a = lParams.weq;
    lLayer.weq_left = @(t, i, W) a(1);
    lLayer.weq_right = @(t, i, W) a(2);
  elseif columns(lParams.weq) == 3
    weq = lParams.weq;
    lLayer.weq_left = @(t, i, W) interp1(weq(:,1), weq(:,2), t);
    lLayer.weq_right = @(t, i, W) interp1(weq(:,1), weq(:,3), t);
  elseif columns(lParams.weq) == 2
    lLayer.weq_left = @(t, i, W) weq(i,1);
    lLayer.weq_right = @(t, i, W) weq(i,2);
  else
    error('Invalid parameter <lParams.weq>');
  endif

  % Matrix geometry contains a struct array
  % .node1 [index]
  % .node2 [index]
  % .r1 [invDiffusion]
  % .r2 [invDiffusion]
  % .geometry1 : geometrical resistance (length over area)
  % .geometry2 :
  lLayer.nodes = arrayfun(@(i) sprintf('n%d',i), 1:lParams.n, 'Uniformoutput', false);
  lLayer.geometry = struct('node1', mat2cell((1:(lParams.n-1))', ones(lParams.n-1,1), 1), ...
                           'node2', mat2cell((2:lParams.n)', ones(lParams.n-1,1), 1),...
                           'invDiffusion1', lLayer.invDiffusion, ...
                           'invDiffusion2', lLayer.invDiffusion, ...
                           'geometry1', mat2cell(lLayer.x_next(1:(end-1)), ones(lParams.n-1,1), 1),
                           'geometry2', mat2cell(lLayer.x_before(2:end), ones(lParams.n-1,1), 1));

  lLayer.command = struct( 'node', {1;  lParams.n}, ...
                           'index', {1; 2},...
                           'invDiffusion',    {lLayer.invDiffusion; lLayer.invDiffusion }, ...
                           'invConductance',  {lLayer.convection_left; lLayer.convection_right }, ...
                           'geometry', { lLayer.x_before(1); lLayer.x_next(end)});

  % Precompute severals vectors to speedup the following computations
  lNode1 = [lLayer.geometry.node1];
  lNode2 = [lLayer.geometry.node2];

  lGeometry1 = [lLayer.geometry.geometry1]';
  lGeometry2 = [lLayer.geometry.geometry2]';
  lCmdNode = [lLayer.command.node];
  lCmdIndex = [lLayer.command.index];
  lCmdGeometry = [lLayer.command.geometry]';

  lNodeVecForMatrix1 = [lNode1 lNode2 lNode1 lNode2 lCmdNode];
  lNodeVecForMatrix2 = [lNode2 lNode1 lNode1 lNode2 lCmdNode];

  Wp = lLayer.winit;
  lR1Vec = lLayer.invDiffusion(lParams.t(1), i, Wp(lNode1));
  lR2Vec = lLayer.invDiffusion(lParams.t(1), i, Wp(lNode2));
  lRNodeCmdVec = lLayer.invDiffusion(lParams.t(1), i, Wp(lCmdNode));
  lRCmdVec = arrayfun(@(v) v.invConductance(lParams.t(1), 1, Wp(v.node)), lLayer.command);
##  In this function, the diffusion function is similar for all nodes
##  so the formulation below is replaced by the formulation above.
##  lR1Vec = arrayfun(@(v) v.invDiffusion1(lParams.t(1), Wp(v.node1)), lLayer.geometry);
##  lR2Vec = arrayfun(@(v) v.invDiffusion2(lParams.t(1), Wp(v.node2)), lLayer.geometry);
##  lR2Vec = lLayer.invDiffusion(lParams.t(1), Wp([lLayer.geometry.node2]));
  lGVec = 1./(lR1Vec .* lGeometry1 + lR2Vec .* lGeometry2);
  lGCmdVec = 1./(lRNodeCmdVec .* lCmdGeometry + lRCmdVec);

  Gp = sparse(lNodeVecForMatrix1, ...
              lNodeVecForMatrix2, ...
              [lGVec; lGVec; -lGVec; -lGVec; -lGCmdVec], ...
              lParams.n, lParams.n);
  Bp = sparse(lCmdNode, ...
              lCmdIndex, ...
              lGCmdVec,
              lParams.n, 2);
  Up = [lLayer.weq_left(lParams.t(1), 1, Wp(1)); lLayer.weq_right(lParams.t(1), 1, Wp(end))];
  lEye = spdiags(lLayer.usize*lParams.length, 0, lParams.n, lParams.n);

  Wmat = NA(lParams.n, numel(lParams.t));
  Wmat(:,1) = lLayer.winit;

  for i=2:numel(lParams.t)
    dt = lParams.t(i)-lParams.t(i-1);

    Wn = Wp; % First pass, coefficients are computed with current water content

    for k=1:2
      % Create the command matrix Bn
      lRCmdVec = arrayfun(@(v) v.invConductance(lParams.t(i), i, Wn(v.node)), lLayer.command);
      lGCmdVec = 1./(lRNodeCmdVec .* lCmdGeometry + lRCmdVec);

      Bn = sparse(lCmdNode, ...
                  lCmdIndex, ...
                  lGCmdVec,
                  lParams.n, 2);

      % Create conductance matrix Gn
      lR1Vec = lLayer.invDiffusion(lParams.t(i), i, Wn(lNode1));
      lR2Vec = lLayer.invDiffusion(lParams.t(i), i, Wn(lNode2));
      lGVec = 1./(lR1Vec .* lGeometry1 + lR2Vec .* lGeometry2);
      Gn = sparse(lNodeVecForMatrix1, ...
                  lNodeVecForMatrix2, ...
                  [lGVec; lGVec; -lGVec; -lGVec; -lGCmdVec], ...
                  lParams.n, lParams.n);

      % Create the command vector
      Un = [lLayer.weq_left(lParams.t(i), i, Wp(1)); lLayer.weq_right(lParams.t(i), i, Wp(end))];

      Ap = lEye+Gp*(dt/2);
      An = lEye-Gn*(dt/2);
      b = Ap*Wp + (Bn*Un + Bp*Up)*dt/2;
      Wtmp = linsolve(An, b, struct('POSDEF', true));
      lMaxChange = max(abs(Wtmp - Wn)./Wn);
      Wn = Wtmp;

      if lOptions.verboseLevel > 1
        disp(sprintf('Iter %03d ; Conv %02d ; Change %.2e', i, k, lMaxChange));
      endif
    endfor

    Wp = Wn;
    Gp = Gn;
    Bp = Bn;

    Wmat(:,i) = Wn;
  endfor

  infos.nodePosition = lLayer.upos;
  infos.nodeVolume = lLayer.usize;

endfunction

