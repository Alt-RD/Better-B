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
% Convert local data to absolute data. The data returned depends on the argument
% <params>.
% Input arguments are:
% 1) F: face object
% 2) <params> = [string] "dims", "center", "corners", "area", "quad"
%   = "dims":    V = [nNodes x 4]
%                Returns the geometry of each nodes of the face in global axis
%   = "center":  V = [3x1]
%                Returns the center of the face
%   = "corners": V = [3x4]  (corners in column u0u0, u1v0, u1v1, u0v1)
%                Returns the position of each 4 corners of the face
%   = "nodesArea":    V = [nNodes x 1]
%                Returns the area of each nodes
%   = "quad":    V = { vertex, face index }
%                  = { [(nxNodes-1)*(nyNodes-1) x 3], [nNodes x 4]
%                Returns geometry information about all nodes of the face. Used
%                to draw nodes in 3D axis.
%   = "resistance":   V = [nNodes x 1]
%                Returns the resistance [K/W] for each node of the face
%   = "rdensity":   V = [nNodes x 1]
%                Returns the resistance density [Km²/W] for each node of the face
%   = "ugrid"/"vgrid": [(n+1) x 1] Returns the grid along the specified direction
%                With n the number of nodes. The grid is the node boundaries.
%   = "xgrid"/"ygrid": [dim 3 x n] Returns a matrix containing the list of 3d
%                vertices along the U (xgrid) or V (ygrid) local coordinate system.
%                The difference with ugrid/vgrid lies in the fact it returns
%                3d vertices in the global axis system.
%   = "xset": [dim n] Returns the set of component values in the u direction in
%                the global space
%   = "vertex": [dim 3 x n] Returns the vertex of a general face in the global axis
function [Ys] = HT_Face_GetAbsoluteData(Fs, params)
  global HT_VAR_EPSILON_POS;

  Ys = cell(numel(Fs),1);

  for i=1:numel(Fs)
    if isstruct(Fs)
      F = Fs(i);
    elseif iscell(Fs)
      F = Fs{i};
    else
      F = Fs;
    endif

    assert(HT_CheckType(F, 'face'), 'Invalid face object');
    lGeneralTypeFace = ~isempty(F.vertex);

##    assert(any(strcmpi(params, {'dims', 'center', 'corners', 'nodesArea', 'quad', 'resistance', ...
##        'ugrid', 'vgrid', 'r', 'xgrid', 'ygrid'})), ...
##      'Function Face_GetAbsoluteData: Invalid parameter specified');
    U=1; V=2;

    if strcmpi(params, 'dims')
      uDir = find(abs(F.axis(:,U)) > 0.5);
      vDir = find(abs(F.axis(:,V)) > 0.5);
      Y = [F.globalPosition(uDir) + (F.axis(uDir, U)*F.size(U))*F.dims(:,[1 2]), ...
           F.globalPosition(vDir) + (F.axis(vDir, V)*F.size(V))*F.dims(:,[3 4])];

      if F.axis(uDir, U) < 0, Y(:,[1 2]) = Y(:,[2 1]); endif;
      if F.axis(vDir, V) < 0, Y(:,[3 4]) = Y(:,[4 3]); endif;
    elseif strcmpi(params, 'center')
      Y = F.globalPosition + F.axis * (0.5*F.size);
    elseif strcmpi(params, 'corners')
      Y = F.globalPosition + F.axis * ([0, 1, 1, 0; 0, 0, 1, 1] .* F.size);
    elseif strcmpi(params, 'nodesArea')
      if ~lGeneralTypeFace % Rectangle based face ?
        Y = (F.dims(:,2) - F.dims(:,1)) .* (F.dims(:,4) - F.dims(:,3)) * prod(F.size);
      else % general face
        %x1 is F.vertex(F.dims(:,1), 1)
        %x2 is F.vertex(F.dims(:,2), 1)
        %y3 is F.vertex(F.dims(:,3), 2)
        % See https://testbook.com/maths/area-of-a-quadrilateral:
        % A=0.5*[(x1y2+x2y3+x3y4+x4y1)-(x2y1+x3y2+x4y3+x1y4)].
        Y = sum(F.vertex(:,1)(F.dims) .* (F.vertex(:,2)(F.dims(:,[2 3 4 1])) - F.vertex(:,2)(F.dims(:,[4 1 2 3]))), 2);
        Y *= 0.5 * prod(F.size);
        assert(all(Y >= 0), 'Invalid order of polygons');
      endif
    elseif strcmpi(params, 'quad')
      error('Not implemented');
    elseif strcmpi(params, 'resistance')
      if isnumeric(F.r)
        lResistance = F.r;
        if numel(lResistance) == 1, lResistance = postpad(lResistance, numel(F.nodes), lResistance); endif;
        Y = lResistance ./ HT_Face_GetAbsoluteData(F, 'area');
      else
        error('Field <r> must be numeric');
      endif
    elseif strcmpi(params, 'r')
      if isnumeric(F.r)
        if isempty(F.r)
          Y = 0;
        else
          Y = F.r;
        endif
##        if numel(Y) == 1, Y = postpad(Y, numel(F.nodes), Y); endif;
      else
        error('Field <r> must be numeric');
      endif
    elseif strcmpi(params, 'ugrid')
      Y = [0; unique(F.dims(:,2))];
      % Make sure there is no numeric error to produce very small nodes
      % In that cases, some fixes will be necessary
      tmp = min(Y(2:end)-Y(1:(end-1)));
      if tmp < 1E-6
        warning(sprintf('Small ugrid step detected <=%e>', tmp));
      endif
    elseif strcmpi(params, 'vgrid')
      Y = [0; unique(F.dims(:,4))];
      tmp = min(Y(2:end)-Y(1:(end-1)));
      if tmp < 1E-6
        warning(sprintf('Small vgrid step detected <=%e>', tmp));
      endif
    elseif strcmpi(params, 'xgrid')
      % xgrid returns a matrix (3 x (nNodes+1)) containing the node boundaries
      % in global coordinate system
      lUGrid = [0; unique(F.dims(:,2))];
      % Make sure there is no numeric error to produce very small nodes
      % In that cases, some fixes will be necessary
      tmp = min(lUGrid(2:end)-lUGrid(1:(end-1)));
      if tmp < 1E-6
        warning(sprintf('Small ugrid step detected <=%e>', tmp));
      endif

      Y = F.globalPosition + (F.axis(:,U) * F.size(U)) .* lUGrid';
    elseif strcmpi(params, 'ygrid')
      % xgrid returns a matrix (3 x (nNodes+1)) containing the node boundaries
      % in global coordinate system
      lVGrid = [0; unique(F.dims(:,4))];
      % Make sure there is no numeric error to produce very small nodes
      % In that cases, some fixes will be necessary
      tmp = min(Y(2:end)-Y(1:(end-1)));
      if tmp < 1E-6
        warning(sprintf('Small ugrid step detected <=%e>', tmp));
      endif

      Y = F.globalPosition + (F.axis(:,V) * F.size(V)) .* lVGrid';
    elseif strcmpi(params, 'xSet')
      if lGeneralTypeFace
##        Y = F.globalPosition + (F.axis .* F.size') * F.vertex';
##        Y = F.axis(:,U)' * Y;
        Y = F.axis(:,U)' * F.globalPosition + F.size(U)*F.vertex(:,1);
        Y = round(Y / HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS;
        Y = unique(Y);
##        Y = uniquetol(Y, HT_VAR_EPSILON_POS, 'DataScale', 1);
      else
        lUGrid = [0; unique(F.dims(:,2))];
##        Y = F.globalPosition + (F.axis(:,U) * F.size(U)) .* lUGrid';
##        Y = F.axis(:,U)' * Y;
        Y = F.axis(:,U)' * F.globalPosition + F.size(U)*lUGrid;
        Y = round(Y / HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS;
        Y = unique(Y);
##        Y = uniquetol(Y, HT_VAR_EPSILON_POS, 'DataScale', 1);
      endif
    elseif strcmpi(params, 'vertex')
      assert(lGeneralTypeFace, 'Not supported for rectangle based face');
      Y = F.globalPosition + (F.axis .* F.size') * F.vertex';
      Y = round(Y / HT_VAR_EPSILON_POS)*HT_VAR_EPSILON_POS;
    else
      error(sprintf('Face_GetAbsoluteData::Invalid parameters <%s>', params));
    endif

    Ys{i} = Y;
  endfor

  if iscell(Fs) || (isstruct(Fs) && (numel(Fs)>1))
  else
    Ys = Ys{1};
  endif
endfunction
