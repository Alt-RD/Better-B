% ========================================================================
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

% @file  HT_Face_Find.m
% @brief Retrieve faces that match a specified criteria
%======================================================================
%
% Find/apply a filter on an array of faces
% Input arguments:
% -> Fs: [structure array] or [cell array] of faces
% -> Supported {parameter, value} pairs
%    "orientation" [dim 3x1] find faces that are oriented in the same direction
%                            The filter relies on <face.norm> vector
%    "returnType" [logical] true -> the function returns a set of faces
%                           false -> the function returns a vector of logical
%    "numberCheck" [scalar] Raise an error if the number of faces that passed
%                           the test is not strictly equal to <numberCheck>
%
% Output:
% if <returnType> is true: a struct array/or cell array containing the 
%                          object that passed the test
% if <returnType> is false: a logical vector specifying which objects
%                           passed the test.
%======================================================================
%> @brief Apply a filter or a test on a set of faces
%>
%> This function finds one or severals faces of a structure array or a cell
%> array.
%> This function supports some input parameters, specified as standard pair
%> ("name", value). Supported parameters are:
%> "orientation" [3xN matrix]: If this parameter is set with a 3d vector D, this 
%>               function returns the faces (or an array of logical values) which are oriented in the same direction.
%>               If a matrix [3xN] is passed, the tested is repeated for each
%>               column and a cell array (dim Nx1) is returned instead. Each cell
%>               contains the set of faces (or logical array) that passed the test.
%> "returnType" [logical]: If false, this function returns a array of logical (size(Y) == size(Fs)).
%>               If true, object faces are returned instead.
%>
%> @param Fs a list of faces (structure array or cell array)
%> @param varargin Parameters passed to the function.
%>
%> @retval Y return value for the first output variable
%======================================================================
function [Y] = HT_Face_Find(Fs, varargin)
  if rows(Fs) == 1, Fs = Fs'; endif
  
  prop = varargin(1:2:end);
  values = varargin(2:2:end);

  assert(iscell(Fs) || isstruct(Fs));
  assert(all(cellfun(@(v) any(strcmpi(v, {'orientation', 'returnType', 'numberCheck'})), prop)), 'Invalid input argument');
  
  lFilter = '';
  lReturnType = true;
  lNumberCheck = NA;
  
  for i=1:numel(prop)
    if strcmpi(prop{i}, 'orientation')
      lFilter = prop{i};
      lVector = values{i};
    elseif strcmpi(prop{i}, 'returnType')
      assert(islogical(values{i}), 'Invalid parameter value');
      lReturnType = values{i};
    elseif strcmpi(prop{i}, 'numberCheck')
      lNumberCheck = values{i};
    else
      error('Invalid parameter');
    endif
  endfor
  
  if strcmpi(lFilter, 'orientation')
    if numel(lVector) == 3
      if iscell(Fs)
        Y = cellfun(@(v) abs(dot(v.norm, lVector)-1) < 1E-10, Fs);
      elseif isstruct(Fs)
        Y = arrayfun(@(v) abs(dot(v.norm, lVector)-1) < 1E-10, Fs);
      endif
    else
      if rows(lVector) != 3, lVector = lVector'; endif
      assert(rows(lVector) == 3, 'Invalid input parmeter');
      
      Y = false(numel(Fs), columns(lVector));
      
      for i=1:columns(lVector)
        if iscell(Fs)
          Y(:,i) = cellfun(@(v) abs(dot(v.norm, lVector(:,i))-1) < 1E-10, Fs);
        elseif isstruct(Fs)
          Y(:,i) = arrayfun(@(v) abs(dot(v.norm, lVector(:,i))-1) < 1E-10, Fs);
        endif
      endfor
    endif
  endif
  
  if ~isna(lNumberCheck)
    assert(all(sum(Y,1) == lNumberCheck), sprintf('Too many object <%d>%d> passed the test', max(sum(Y,1)), lNumberCheck));
  endif
  
  if lReturnType
    if iscolumn(Y) == 1
      Y = Fs(Y);
    else
      lBoolMatrix = Y;
      Y = cell(1, columns(lBoolMatrix));
      
      for i=1:columns(lBoolMatrix)
        Y{1,i} = Fs(lBoolMatrix(:,i));
      endfor
    endif
  endif
  
endfunction
