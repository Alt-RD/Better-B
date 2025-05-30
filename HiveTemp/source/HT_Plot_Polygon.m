%  This file is part of project HiveTemp.
%
%  Copyright (c) 2023 AltRD-Emmanuel Ruffio
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
function [V F h] = HT_Plot_Polygon(polygonType, varargin)
  props = varargin(1:2:end);
  values = varargin(2:2:end);

  params = struct();
  for i=1:numel(props)
    params = setfield(params, props{i}, values{i});
  endfor

  if strcmpi(polygonType, 'disk')
    params = HT_CheckField(params, 'n',            32);
    params = HT_CheckField(params, 'radius',       1);
    params = HT_CheckField(params, 'position',     zeros(3,1));
    params = HT_CheckField(params, 'draw',         false);
    params = HT_CheckField(params, 'color',        [1 1 1]);
    params = HT_CheckField(params, 'edgeColor',    [0 0 0]);

    if isscalar(params.n)
      t = (0:params.n)'/params.n*2*pi;
      ny = int32(params.n);
    else
      t = params.n;
      assert(all(diff(t) > 0), 'Invalid field <n>. Must be a scalar or a vector V with V(end)==V(1)');
      ny = numel(t)-1;
    endif

    x = params.radius(1)*cos(t);
    y = params.radius(end)*sin(t);

    V = [x, y, zeros(size(x))];
    V = [0, 0, 0; V];
    V += params.position(:)';

    F = [ones(ny,1), (2:(ny+1))', (3:(ny+2))'];
    h = 0;

    if params.draw
      h = patch("Faces", F, "Vertices", V, ...
                'edgecolor', params.edgeColor, ...
                'facealpha', 1.0, ...
                'facelighting', 'none', ...
                'facecolor', params.color, ...
                'faceVertexCData', [1 0 0]);
    endif
  elseif strcmpi(polygonType, 'cylinder')
    params = HT_CheckField(params, 'n',            32);
    params = HT_CheckField(params, 'radius',       1);
    params = HT_CheckField(params, 'height',       1);
    params = HT_CheckField(params, 'position',     zeros(3,1));
    params = HT_CheckField(params, 'draw',         false);
    params = HT_CheckField(params, 'color',        [1 1 1]);
    params = HT_CheckField(params, 'edgeColor',    [0 0 0]);
    params = HT_CheckField(params, 'full',         true); % Add top and bottom face
    params = HT_CheckField(params, 'order',        'ccw', @(v) any(strcmpi(v, {'ccw', 'cw'}))); % Counter clock wise
    params = HT_CheckField(params, 'sort',         'theta', @(v) any(strcmpi(v, {'theta', 'height'})));

    if isscalar(params.height), params.height = [0; params.height]; endif;
    if isscalar(params.radius), params.radius = repmat(params.radius, numel(params.height), 1); endif
    if ~iscolumn(params.height) params.height = params.height'; endif;
    if ~iscolumn(params.radius) params.radius = params.radius'; endif;
    %if isrow(params.color) params.color = params.color'; endif;
    if iscolumn(params.color) params.color = repmat(params.color, 1, 3); endif

    assert(numel(params.radius) == numel(params.height));
    assert(iscolumn(params.height));
    assert(iscolumn(params.radius));
    assert(strcmpi(params.order, 'ccw'), 'Not implemented. <order> must be ccw');
    assert(columns(params.color) == 3, 'Invalid field <color>. Matrix has not a valid size');
##    assert(~params.full, 'Not implemented');

##    assert(params.draw, 'Not supported yet. draw must be true');
##
##    [x,y,z] = cylinder([params.radius(1) params.radius(end)], params.n); hold on;
##    z *= params.height;
##    x += params.position(1);
##    y += params.position(2);
##    z += params.position(3);
##    hm = surf(x,y,z);
##    set(hm, 'edgecolor', params.edgeColor, ...
##            'facecolor', params.color);

    if isscalar(params.n)
      t = (0:(params.n-1))'/params.n*2*pi;
      ny = int32(params.n);
    else
      t = params.n;
      if abs(mod(t(end),2*pi)-mod(t(1),2*pi)) < eps('double'), t(end) = []; endif

      assert(all(diff(t) > 0), 'Invalid field <n>. Must be a scalar or a vector V with V(end)==V(1)');
      ny = numel(t);
    endif

    nz = int32(numel(params.height)-1);

    x = cos(t);
    y = sin(t);

    V = [ (x * params.radius')(:), (y * params.radius')(:), repelem(params.height, ny, 1)];

    if strcmpi(params.sort, 'height')
      F = [(0:(nz-1))'*ny, 1+(0:(nz-1))'*ny, 1+(1:nz)'*ny, (1:nz)'*ny];
      F = repmat(F, ny, 1) + repelem((0:(ny-1))', nz, 1);
      F((end-nz+1:end) ,2:3) -= ny; % Vertices of the last vertical (theta=2pi) are in fact the first vertical (theta=0)
    else
      F = [(0:(ny-1))' , mod((1:ny)', ny)];
      F = [F, fliplr(F)+ny];

      if nz > 1
        F = repmat(F, nz, 1) + repelem(ny*(0:(nz-1))' , ny, 1);
      endif
    endif

    if params.full
      V = [V; 0, 0, params.height(1); 0, 0, params.height(end)];
      iZeroTop = rows(V)-1;
      iZeroBot = rows(V)-2;

      Fadd = ny*nz + [(0:(ny-1))', [(1:(ny-1)), 0]'];
      Fadd = [Fadd, repmat(iZeroTop, ny, 2)];
      F = [F; Fadd];

      Fadd = [(0:(ny-1))', [(1:(ny-1)), 0]'];
      Fadd = [Fadd, repmat(iZeroBot, ny, 2)];
      F = [F; Fadd];
    endif

    V += params.position(:)';
    F += 1;

    assert(any(rows(params.color) == [1 rows(F)]), sprintf('Invalid color matrix <color>. Must be 1x3 or Nfaces(=%d)x3', rows(F)));

    if params.draw
      lFaceVertexCData = [];
      lFaceColor = params.color;

      if ~isrow(params.color)
        lFaceColor = 'flat';
        lFaceVertexCData = params.color;
      endif

      h = patch("Faces", F, "Vertices", V, ...
                'edgecolor', params.edgeColor, ...
                'facealpha', 1.0, ...
                'facelighting', 'none', ...
                'facecolor', lFaceColor, ...
                'facevertexcdata', lFaceVertexCData);
    endif
  else
    error('Invalid polygon type');
  endif
endfunction
