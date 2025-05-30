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
%> @class HT_Rect2
%> @brief Defines tools to manage 2D rectangles
%>
%======================================================================
classdef HT_Rect2
  properties
    m_origin = zeros(2,1);
    m_axis = [1 0; 0 1]';
    m_size = [];
  endproperties

  methods (Access = public)
    function R = HT_Rect2(_origin, _axis, _size)
      if nargin == 0
        return;
      endif

      assert(nargin == 3);
      assert(numel(_size) == 2);

      % Create a 2D rect from 2D data
      if numel(_origin) == 2
        assert(all(size(_axis) == [2 2]));

        R.m_origin = _origin(:);
        R.m_axis = _axis ./ vecnorm(_axis);
        R.m_size = _size(:);
      else
        error('Invalid input argument');
      endif
    endfunction

    function disp(obj)
      disp(sprintf('HT_Rect2: o=[%#+.3G;%#+.3G], u=[%#+.2f;%#+.2f], v=[%#+.2f;%#+.2f], size=[%#+.3G;%#+.3G]', obj.m_origin, obj.m_axis, obj.m_size));
    endfunction

    function B = isNull(obj)
      B = isempty(obj.m_size);
    endfunction

    function B = isEmpty(obj)
      B = (prod(obj.m_size) < 1E-10);
    endfunction

    function V = size(obj)
      V = obj.m_size;
    endfunction

    function V = origin(obj)
      V = obj.m_origin;
    endfunction

    function V = axis(obj)
      V = obj.m_axis;
    endfunction

    function R = adaptAxis(obj, arg)
      if strcmp(class(arg), 'HT_Rect2')
        arg = arg.m_axis;
      endif

      assert(all(size(arg) == [ 2 2]));
      M = obj.m_axis' * arg;
      assert(all(any(abs(M(:) - [-1, 0, 1]) < 1E-10, 2)), 'Rect2 objects have not compatible axis');

      if abs(M(1,1)) < 1E-10
        obj.m_size = flip(obj.m_size);
        obj.m_axis = fliplr(obj.m_axis);

        M = flipud(M);
      endif

      % If uA and uB are reversed, uB is reversed
      if M(1,1) < 0
        obj.m_origin += obj.m_size(1) * obj.m_axis(:,1);
        obj.m_axis(:,1) = - obj.m_axis(:,1);
      endif

      % If vA and vB are reversed, vB is reversed
      if M(2,2) < 0
        obj.m_origin += obj.m_size(2) * obj.m_axis(:,2);
        obj.m_axis(:,2) = - obj.m_axis(:,2);
      endif

      R = obj;
    endfunction

    %
    %======================================================================
    %> @brief Returns the intersection of rects
    %======================================================================
    function R = intersectRect(obj, argList)
      if ~iscell(argList)
        argList = { argList };
      endif

      assert(~isempty(argList), 'Rect2::intersectRect: Missing arguments');

      R = HT_Rect2(); % Default is null object

      for i=1:numel(argList)
        arg = argList{i};
        assert(strcmp(class(arg), 'HT_Rect2'));

        arg = arg.adaptAxis(obj);

        % Shrinks the rect <arg> to adapt to <obj>
        lArgMove = arg.m_axis' * (obj.m_origin - arg.m_origin);
        % Only positive move is allowed, it means the rect can only be reduced
        lArgMove = min(max(lArgMove, 0), arg.m_size);
        arg.m_size -= lArgMove;
        arg.m_origin += arg.m_axis*lArgMove;

        % The same operation is applied on rect 1 <obj>
        lObjMove = arg.m_axis' * (arg.m_origin - obj.m_origin);
        lObjMove = min(max(lObjMove, 0), obj.m_size);
        obj.m_size -= lObjMove;
        obj.m_origin += obj.m_axis*lObjMove;

        % If the two origin could not be put in contact, it means there is
        % no intersection.
        if norm(obj.m_origin - arg.m_origin, Inf) > 1E-10 % No intersection
          return;
        endif

        lObjOppositeOrigin = obj.oppositeOrigin();

        lArgMove = arg.m_axis' * (lObjOppositeOrigin - arg.oppositeOrigin());
        lArgMove = max(min(lArgMove, 0), -arg.m_size);
        arg.m_size += lArgMove;

        % The same operation is applied on rect 1 <obj>
        lObjMove = arg.m_axis' * (arg.oppositeOrigin() - lObjOppositeOrigin);
        lObjMove = max(min(lObjMove, 0), -obj.m_size);
        obj.m_size += lObjMove;
      endfor

      R = obj;
    endfunction

    %
    %======================================================================
    %> @brief Returns the union of rects
    %> @image html Rect2_UnionRect_Examples.png
    %!["Example"](./../images/Rect2_UnionRect_Examples.png)
    %======================================================================
    function R = unionRect(obj, argList, varargin)
      if ~iscell(argList)
        argList = { argList };
      endif

      % This operation is usefull when re-calling this function with a previous
      % one that takes varargin argument as input. It allows to call this function
      % like unionRect(myObj, myList, varargin).
      if (nargin == 3) && iscell(varargin) && (numel(varargin) == 1)
        varargin = varargin{1};
      endif

      R = HT_Rect2(); % Default is null object

      lAllowIntersect = true;
      lAllowEmptySpace = true;
      lNoError = false;

      prop = varargin(1:2:end);
      value = varargin(2:2:end);

      assert(numel(prop) == numel(value), 'Invalid input parameters');

      for i=1:numel(prop)
        if strcmpi(prop{i}, 'allowEmptySpace')
          lAllowEmptySpace = value{i};
        elseif strcmpi(prop{i}, 'allowIntersect')
          lAllowIntersect = value{i};
        elseif strcmpi(prop{i}, 'noerror')
          lNoError = true;
        else
          error(sprintf('Invalid input parameter <%s>', prop{i}));
        endif
      endfor

      for i=1:numel(argList)
        arg = argList{i};
        assert(strcmp(class(arg), 'HT_Rect2'));

        arg = arg.adaptAxis(obj);

        lBoundingOrigin = min(obj.m_origin, arg.m_origin);
        lBoundingOppOrigin = max(obj.oppositeOrigin(), arg.oppositeOrigin());

        lBounding = HT_Rect2( lBoundingOrigin, ...
                              obj.m_axis, ...
                              lBoundingOppOrigin - lBoundingOrigin);

        if ~lAllowEmptySpace
          if ~all(lBounding.onBoundaries(obj.corners()))
            assert(lNoError, 'Rect2::Union: Rects are not aligned');
            return;
          endif
          if ~all(lBounding.onBoundaries(arg.corners()))
            assert(lNoError, 'Rect2::Union: Rects are not aligned');
            return;
          endif

          lBArea = lBounding.area();
          lObjArea = obj.area();
          lArgArea = arg.area();

          if (lBArea > lObjArea + lArgArea + 1E-10)
            assert(lNoError, 'Rect2::Union: Rects are disjoint');
            return;
          endif

          if ~lAllowIntersect && (lBArea < lObjArea + lArgArea - 1E-10)
            assert(lNoError, 'Rect2::Union: Rects are intersecting');
          endif
        else
          if ~lAllowIntersect
            if obj.intersectRect(arg).area() > 1E-10
              assert(lNoError, 'Rect2::Union: Intersection is not empty');
              return;
            endif
          endif
        endif

        obj = lBounding;
      endfor

        R = obj;
    endfunction

    %> @brief Returns the position of a list of points in the local axis of the object
    %> @param _points (dim 2xN)
    function P = toLocalAxis(obj, _points)
      assert(rows(_points) == 2);
      P = obj.m_axis' * (_points - obj.m_origin);
    endfunction

    function C = corners(obj)
      C = obj.m_origin + obj.m_axis * ([0 0 1 1;0 1 1 0] .* obj.m_size);
    endfunction

    function B = contains(obj, _points)
      _points = obj.toLocalAxis(_points);
      B = (_points(1,:) <= obj.m_size(1) + 1E-10) & (-1E-10 <= _points(1,:)) & ...
          (_points(2,:) <= obj.m_size(2) + 1E-10) & (-1E-10 <= _points(2,:));
    endfunction

    function B = onBoundaries(obj, _points)
      _points = obj.toLocalAxis(_points);
      B = any(abs([0; obj.m_size(1)] - _points(1,:)) < 1E-10) | any(abs([0; obj.m_size(2)] - _points(2,:)) < 1E-10);
    endfunction

    function P = oppositeOrigin(obj)
      P = obj.m_origin + obj.m_axis * obj.m_size;
    endfunction

    function V = area(obj)
      V = prod(obj.m_size);
    endfunction

    function V = to3dSpace(obj, _points, _axis, _origin)
      if nargin < 4, _origin = zeros(3,1); endif;

      V = dot(_axis(:,3),_origin)*_axis(:,3) + _axis(:,1:2) * _points;
    endfunction

    function plot(obj, varargin)
      prop = varargin(1:2:end);
      value = varargin(2:2:end);

      lIsPatch = false;

      assert(numel(prop) == numel(value), 'Invalid parameters');
      for i=1:numel(prop)
        if strcmpi(prop{i}, 'faceColor')
          lIsPatch = true;
          break;
        endif
      endfor

      rectangle('Position', [obj.m_origin'-0.05*mean(obj.m_size), repmat(0.1*mean(obj.m_size), 1, 2) ], 'FaceColor', 'black');
      %h = rectangle('Position', [obj.m_origin', obj.m_size']);

      C = obj.corners();

      if lIsPatch
        h = patch(C(1,:), C(2,:));
      else
        C = [C, C(:,1)];
        h = line(C(1,:), C(2,:));
      endif

      assert(numel(prop) == numel(value), 'Invalid parameters');
      for i=1:numel(prop)
        set(h, prop{i}, value{i});
      endfor
    endfunction

    function plot3d(obj, _axis, _origin, varargin)
      if nargin < 3, _origin = zeros(3,1); endif;

      prop = varargin(1:2:end);
      value = varargin(2:2:end);

      lIsPatch = false;

      assert(numel(prop) == numel(value), 'Invalid parameters');
      for i=1:numel(prop)
        if strcmpi(prop{i}, 'faceColor')
          lIsPatch = true;
          break;
        endif
      endfor

      C = dot(_axis(:,3),_origin)*_axis(:,3) + _axis(:,1:2) * obj.corners();

      if lIsPatch
        h = patch(C(1,:), C(2,:), C(3,:));
      else
        C = [C, C(:,1)];
        h = line(C(1,:), C(2,:), C(3,:));
      endif

      assert(numel(prop) == numel(value), 'Invalid parameters');
      for i=1:numel(prop)
        set(h, prop{i}, value{i});
      endfor
    endfunction
  endmethods

  methods (Static = true)
    function R = from3dSpace(_origin, _axis, _size, _system, _systemOrigin)
      if nargin < 5, _systemOrigin = zeros(3,1); endif;

      assert(numel(_systemOrigin) == 3);
      assert(numel(_origin) == 3);
      assert(all(size(_axis) == [3 2]));
      assert(numel(_size) == 2);
      assert(norm(_axis' * _system(:,3), Inf) < 1E-10, 'Rect2::from3dSpace: specified axis do not belong to plane xy');

      _origin = _origin(:);
      _systemOrigin = _systemOrigin(:);

      % Change coordinate a and remove z component
      R = HT_Rect2((_system' * (_origin-_systemOrigin))(1:2), (_system' * _axis)(1:2, :), _size);
    endfunction

    function R = fromUnion(_list, varargin)
      assert(iscell(_list));

      if isempty(_list)
        R = HT_Rect2();
      else
        R = _list{1}.unionRect(_list(2:end), varargin);
     endif
    endfunction

    function R = fromIntersection(_list)
      assert(iscell(_list));

      if isempty(_list)
        R = HT_Rect2();
      else
        R = _list{1}.intersectRect(_list(2:end));
     endif
    endfunction
  endmethods
endclassdef
