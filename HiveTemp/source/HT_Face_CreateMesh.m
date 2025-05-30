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

%======================================================================
%> @brief Brief Returns a copy of the face with a mesh
%>
%> @param face an object of type <face>
%>
%> @retval F a new face
%======================================================================
function Fs = HT_Face_CreateMesh(Fs, varargin)
  global HT_VAR_EPSILON_POS;
  global HT_VAR_EPSILON_U;

  assert(nargin >= 1, "Missing argument");
  assert(mod(numel(varargin), 2) == 0, "Invalid input parameters");

  prop = varargin(1:2:end);
  value = varargin(2:2:end);

  lGridType = 'uniform';
  lUEdgeType = 'ff';
  lVEdgeType = 'ff';
  lN = [];

  for i=1:numel(prop)
    if strcmpi(prop{i}, 'n')
      assert(isnumeric(value{i}) && any(numel(value{i}) == [1 2]) && all((value{i} > 0) && (value{i} <= 1000)), 'Invalid parameter <n>');
      lN = int32(value{i});
    elseif strcmpi(prop{i}, 'type')
      assert(any(strcmpi(value{i}, 'uniform')), sprintf('Invalid mesh type <%s>', value{i}));
      lGridType = value{i};
    elseif strcmpi(prop{i}, 'uEdgeType')
      assert(any(strcmpi(value{i}, {'ff', 'hf', 'fh', 'hh'})), sprintf('Invalid u edge type <%s>', value{i}));
      lUEdgeType = value{i};
    elseif strcmpi(prop{i}, 'vEdgeType')
      assert(any(strcmpi(value{i}, {'ff', 'hf', 'fh', 'hh'})), sprintf('Invalid v edge type <%s>', value{i}));
      lVEdgeType = value{i};
    else
      error(sprintf('Invalid input parameter name <%s>', prop{i}));
    endif
  endfor

  assert(all(lN > 0), 'Invalid node count');
  assert(~isempty(lN), 'Node count, parameter <n> is not defined');
  lN = postpad(lN, 2, lN(1));

  du = 1 / (double(lN(1)) - 0.5*(lUEdgeType(1)=='h') - 0.5*(lUEdgeType(2)=='h'));
  dv = 1 / (double(lN(1)) - 0.5*(lVEdgeType(1)=='h') - 0.5*(lVEdgeType(2)=='h'));

  lUSizeVec = repmat(du, lN(1), 1);
  lVSizeVec = repmat(dv, lN(2), 1);

  if (lUEdgeType(1)=='h'), lUSizeVec(1) *= 0.5; endif
  if (lUEdgeType(2)=='h'), lUSizeVec(end) *= 0.5; endif
  if (lVEdgeType(1)=='h'), lVSizeVec(1) *= 0.5; endif
  if (lVEdgeType(2)=='h'), lVSizeVec(end) *= 0.5; endif

  lUDim = [0; cumsum(lUSizeVec)];
  lVDim = [0; cumsum(lVSizeVec)];

  lUPos = 0.5*(lUDim(1:(end-1)) + lUDim(2:end));
  if (lUEdgeType(1)=='h'), lUPos(1) = 0; endif
  if (lUEdgeType(2)=='h'), lUPos(end) = 1; endif

  lVPos = 0.5*(lVDim(1:(end-1)) + lVDim(2:end));
  if (lVEdgeType(1)=='h'), lVPos(1) = 0; endif
  if (lVEdgeType(2)=='h'), lVPos(end) = 1; endif

  lUPos = round(lUPos / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  lVPos = round(lVPos / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  lUDim = round(lUDim / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  lVDim = round(lVDim / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  lUSizeVec = round(lUSizeVec / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;
  lVSizeVec = round(lVSizeVec / HT_VAR_EPSILON_U) * HT_VAR_EPSILON_U;

  lPos = [repmat(lUPos, lN(2), 1), repelem(lVPos, lN(1))];
  lDims = [repmat([lUDim(1:(end-1)), lUDim(2:end)], lN(2), 1), ...
            repelem(lVDim(1:(end-1)), lN(1), 1), ...
            repelem(lVDim(2:end), lN(1), 1)];

  lNodes = cell(lN(1)*lN(2), 1);
    for m=1:lN(2)
      lNodes((1+(m-1)*lN(1):(m*lN(1)))) = ...
        arrayfun(@(v) sprintf('.n%d_%d', v, m), [1:lN(1)], 'UniformOutput', false);
    endfor

  for i=1:numel(Fs)
    F = Int_GetObject(Fs(i));

    assert(~isempty(F.size) && all(F.size > 0), 'Face size is not defined');

    F.n = lN;
    F.pos = lPos;
    F.dims = lDims;
    F.nodes = cellfun(@(v) strcat(F.name, v), lNodes, 'UniformOutput', false);

    if iscell(Fs)
      Fs{i} = F;
    elseif isstruct(Fs)
      Fs(i) = F;
    else
      error('Invalid face object');
    endif
  endfor

endfunction

% Return an object. If obj is a structure, it does nothing.
% If a cell, it returns the content.
function obj = Int_GetObject(obj)
  if iscell(obj)
    obj = obj{1};
  endif
endfunction
