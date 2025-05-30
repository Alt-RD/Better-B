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
%> @class HT_Quad2
%> @brief Defines tools to manage 2D quadrilateral
%>
%======================================================================
classdef HT_Quad2
  properties
    m_vx = [];
    m_vy = [];
    m_bd = [];
  endproperties

  methods (Access = public)
    function R = HT_Quad2(varargin)
      if nargin == 0
        return;
      elseif nargin == 1
        _vxy = varargin{1};
        if columns(_vxy) != 4, _vxy = _vxy'; endif;

        assert(all(size(_vxy) == [2 4]), 'HT_Quad2: Invalid argument <x>');
        R.m_vx = _vxy(1,:)';
        R.m_vy = _vxy(2,:)';
      elseif nargin == 2
        _vx = varargin{1}
        _vy = varargin{2}

        assert(isnumeric(_vx), 'HT_Quad2: Invalid argument <x>');
        assert(isnumeric(_vy), 'HT_Quad2: Invalid argument <y>');

        if isscalar(_vx)
          R.m_vx = _vx * ones(4,1);
        elseif numel(_vx) == 3
          assert(numel(_vx) == 4, 'HT_Quad2: Invalid size of argument <x>');
          R.m_vx = _vx(:);
        endif

        if isscalar(_vy)
          R.m_vy = _vy * ones(4,1);
        elseif numel(_vy) == 3
          assert(numel(_vy) == 4, 'HT_Quad2: Invalid size of argument <y>');
          R.m_vy = _vy(:);
        endif
      else
        error('HT_Quad2: Invalid input argument');
      endif

      m_bd = boundingBox([R.m_vx', R.m_vy']);
    endfunction

    function disp(obj)
      if obj.isNull()
        disp('HT_Quad2: null');
      else
        disp(sprintf('HT_Quad2: area=%.3G, v1=[%#+.3G;%#+.3G], v2=[%#+.3G;%#+.3G], v3=[%#+.3G;%#+.3G], v4=[%#+.3G;%#+.3G]', ...
             obj.area(), [obj.m_vx obj.m_vy]'(:)));
      endif
    endfunction

    function B = isNull(obj)
      B = isempty(obj.m_vx) && isempty(obj.m_vy);
    endfunction

    function B = isEmpty(obj)
      B = (obj.area() < 1E-10);
    endfunction

    % See https://testbook.com/maths/area-of-a-quadrilateral:
    % A=0.5*[(x1y2+x2y3+x3y4+x4y1)-(x2y1+x3y2+x4y3+x1y4)].
    function V = area(obj)
      if (numel(obj.m_vx) != 4) || (numel(obj.m_vy) != 4)
        V = 0;
      else
        V = 0.5 * sum(obj.m_vx .* (obj.m_vy([2 3 4 1]) - obj.m_vy([4 1 2 3])));
      endif
    endfunction

    % Returns the corners in a matrix [2x4] = [u1 u2 u3 u4; v1 v2 v3 v4];
    function C = corners(obj)
      C = [obj.m_vx'; obj.m_vy'];
    endfunction

    % Returns the vector associated to each segment
    function C = vectors(obj, dim)
      if nargin < 2, dim = 2; endif;

      C = [ obj.m_vx([2 3 4 1]) - obj.m_vx([1 2 3 4]) , obj.m_vy([2 3 4 1]) - obj.m_vy([1 2 3 4]) ]';
      if (dim == 3), C = [C; zeros(1,4)]; endif;
    endfunction

    % Returns the vector associated to each segment
    function C = normalizedVectors(obj, dim)
      if nargin < 2, dim = 2; endif;

      C = [ obj.m_vx([2 3 4 1]) - obj.m_vx([1 2 3 4]) , obj.m_vy([2 3 4 1]) - obj.m_vy([1 2 3 4]) ]';
      if (dim == 3), C = [C; zeros(1,4)]; endif;
      C ./= vecnorm(C, 2, 1);
    endfunction

    % Returns the normals to each segment, towards the quad interior
    % if vertices are defined using counter clockwise order;
    % Direction forces the normal to point outside (direction == 1)
    % or inside (direction = -1), or related to the vertex order (direction = [] or 0)
    function C = normals(obj, direction, dim)
      if nargin < 2, direction = []; endif % do not check direcition by default
      if nargin < 3, dim = 2; endif

      V = obj.normalizedVectors(3);
      C = crossProduct3d([0,0,1], V')';

      if ~isempty(direction) && any(direction == [-1 1])
        lOrder = obj.order(V);
        assert(lOrder != 0, 'HT_Quad2: unable to retrieve quad order');
        if (direction != lOrder)
          C = -C;
        endif
      endif

      if (dim == 2), C(3,:) = []; endif
    endfunction

    % Returns the order used to define the quad
    % 1: clockwise / 0: undefined / -1:counter clockwise
    function C = order(obj, v)
      if nargin < 2
        V = obj.vectors(3);
      else
        V = v;
      endif
      V = [V, V(:,1)];
      C = cross(V(:,1:(end-1)), V(:,2:end));
      C = C(3,:);

      if all(abs(C) < 1E-10)
        C = 0;
      elseif all(C >= -1E-10)
        C = -1;
      elseif all(C <= 1E-10)
        C = 1;
      else
        C = 0;
      endif
    endfunction

    %======================================================================
    %> @brief Test wether or not a quad contains a set of points
    %======================================================================
    function B = containsPoints(obj, _points)
      if (nargin < 2) || isempty(_points)
        B = [];
        return;
      endif

      if columns(_points) != 2, _points = _points'; endif;
      assert(columns(_points) == 2, 'HT_Quad: Invalid point matrix');

      nPoint = rows(_points);
      C = repelem(_points, 4,1) - repmat([obj.m_vx , obj.m_vy], nPoint, 1);
      N = obj.normals(-1, 2);

      C = dot(C, repmat(N', nPoint, 1), 2);
      C = reshape(C, 4, nPoint);
      B = all(C >= -1E-10, 1)';
    endfunction

    %======================================================================
    %> @brief Test wether or not a quad contains a set of quads
    %======================================================================
    function B = containsQuads(obj, _quadList)
      B = cellfun(@(q) all(obj.containsPoints(q.corners())), _quadList);
    endfunction
##      _points = obj.toLocalAxis(_points);
##      B = (_points(1,:) <= obj.m_size(1) + 1E-10) & (-1E-10 <= _points(1,:)) & ...
##          (_points(2,:) <= obj.m_size(2) + 1E-10) & (-1E-10 <= _points(2,:));

    function B = onBoundaries(obj, _points)
##      _points = obj.toLocalAxis(_points);
##      B = any(abs([0; obj.m_size(1)] - _points(1,:)) < 1E-10) | any(abs([0; obj.m_size(2)] - _points(2,:)) < 1E-10);
    endfunction

##    function plot(obj, varargin)
##      prop = varargin(1:2:end);
##      value = varargin(2:2:end);
##
##      lIsPatch = false;
##
##      assert(numel(prop) == numel(value), 'Invalid parameters');
##      for i=1:numel(prop)
##        if strcmpi(prop{i}, 'faceColor')
##          lIsPatch = true;
##          break;
##        endif
##      endfor
##
##      rectangle('Position', [obj.m_origin'-0.05*mean(obj.m_size), repmat(0.1*mean(obj.m_size), 1, 2) ], 'FaceColor', 'black');
##      %h = rectangle('Position', [obj.m_origin', obj.m_size']);
##
##      C = obj.corners();
##
##      if lIsPatch
##        h = patch(C(1,:), C(2,:));
##      else
##        C = [C, C(:,1)];
##        h = line(C(1,:), C(2,:));
##      endif
##
##      assert(numel(prop) == numel(value), 'Invalid parameters');
##      for i=1:numel(prop)
##        set(h, prop{i}, value{i});
##      endfor
##    endfunction

##    function plot3d(obj, _axis, _origin, varargin)
##      if nargin < 3, _origin = zeros(3,1); endif;
##
##      prop = varargin(1:2:end);
##      value = varargin(2:2:end);
##
##      lIsPatch = false;
##
##      assert(numel(prop) == numel(value), 'Invalid parameters');
##      for i=1:numel(prop)
##        if strcmpi(prop{i}, 'faceColor')
##          lIsPatch = true;
##          break;
##        endif
##      endfor
##
##      C = dot(_axis(:,3),_origin)*_axis(:,3) + _axis(:,1:2) * obj.corners();
##
##      if lIsPatch
##        h = patch(C(1,:), C(2,:), C(3,:));
##      else
##        C = [C, C(:,1)];
##        h = line(C(1,:), C(2,:), C(3,:));
##      endif
##
##      assert(numel(prop) == numel(value), 'Invalid parameters');
##      for i=1:numel(prop)
##        set(h, prop{i}, value{i});
##      endfor
##    endfunction
    function h = draw(obj)
      h = drawPolygon(obj.m_vx, obj.m_vy);
    endfunction
  endmethods

  methods (Static = true)
    % Convert quads in 3d space into Quad2 object based on plane definition
    % (_axis and _systemOrigin)
    % _vertex is [3 x 4N] with N the number of quads to generate
    % Returns a cell array of Quad2 objects
    function R = from3dSpace(_vertex, _axis, _systemOrigin)
      if nargin < 3, _systemOrigin = zeros(3,1); endif;

      assert(rows(_vertex) == 3);
      assert(mod(columns(_vertex),4) == 0);
      assert(all(size(_systemOrigin) == [3 1]));
      assert(all(size(_axis) == [3 2]));

      nQuads = columns(_vertex) / 4;
      V = _axis' * (_vertex - _systemOrigin);

      R = cell(nQuads, 1);
      R = arrayfun(@(i) HT_Quad2(V(:,4*(i-1)+1+(0:3))) , 1:nQuads, 'uniformOutput', false);
    endfunction

##    function R = fromUnion(_list, varargin)
##      assert(iscell(_list));
##
##      if isempty(_list)
##        R = HT_Rect2();
##      else
##        R = _list{1}.unionRect(_list(2:end), varargin);
##     endif
##    endfunction
##
##    function R = fromIntersection(_list)
##      assert(iscell(_list));
##
##      if isempty(_list)
##        R = HT_Rect2();
##      else
##        R = _list{1}.intersectRect(_list(2:end));
##     endif
##    endfunction
  endmethods
endclassdef
