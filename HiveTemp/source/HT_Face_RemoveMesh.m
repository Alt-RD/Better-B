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
%> @brief Brief Returns a copy of the face with no mesh
%>
%> This function removes all information related to nodes and meshes.
%> The new face has only 1 nodes at the center of the geometry
%> Information about thermal resistance is removed as well.
%>
%> @param face an object of type <face>
%>
%> @retval F a new face
%======================================================================
function F = HT_Face_RemoveMesh(face)
  assert(HT_CheckType(face, 'face'), 'Invalid face object');

  for i=1:numel(face)
    f = Int_GetObject(face,i);
    f = HT_Face_Init(f.name,  'copy', f,             ...
                              'dims', [],     ...
                              'nodes', {},        ...
                              'position', [],    ...
                              'material', [], ...
                              'materialIndex', [], ...
                              'n', [],               ...
                              'r', []);
##    f = HT_Face_Init(f.name,  'copy', f,             ...
##                              'dims', [0, 1, 0, 1],     ...
##                              'nodes', {'node'},        ...
##                              'position', [0.5 0.5],    ...
##                              'n', [1 1],               ...
##                              'r', []);

    if iscell(face)
      F{i} = f;
    else
      F(i) = f;
    endif
  endfor

endfunction

function M = Int_GetObject(X, ind)
  if nargin < 2, ind = 1; endif;

  if iscell(X)
    M=X{ind};
  else
    M=X(ind);
  endif
endfunction

