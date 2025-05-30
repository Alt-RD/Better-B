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
% This function returns a structure containing information about the shape
% used to defines a face.
% It is used mainly to speed up face intersections and to do other stuff like this.
%
% Most field are builtin, but fields can be added freely as long as they don't
% interfere with buitlin fields
% Details:
% Shape: "circle":
% -> uSize : [Nx1] with N > 1 for multilayer circle and/or hollow circle
% -> vSize : [Nx1]
% -> uGrid: [(N+1)x1] with N the number of nodes along the radius axis
% -> vGrid: [Nx1] with N the number the number of subdivision of the circle
%           [1x1 integer] the number of subdivision of the circle
% -> uNode: [Nx1] with N the number of layer. Each component is the number of
%                 nodes in that layer
% -> uNodePos: [Nx1] the relative position of each nodes (relative to the thickness)
%                    0 in the internal radius and 1 the outer radius
function M = HT_Face_InitMetaData(varargin)
  props = varargin(1:2:end);
  values = varargin(2:2:end);

  assert(numel(props) == numel(values), 'Invalid argument');
  lTypeIndex = find(strcmpi(props, 'type'));
  assert(~isempty(lTypeIndex), 'Type shape must be defined');
  assert(any(strcmpi(values{lTypeIndex}, 'circle')), sprintf('Unknown shape <%s>', values{lTypeIndex}));

  M = struct( 'type', values{lTypeIndex}, ...
              'uSize', [], ...
              'vSize', [], ...
              'uGrid', [], ...
              'vGrid', [], ...
              'uNode', [], ...
              'uNodePos', []);

  if strcmpi(M.type, 'circle')
    M = Int_LoadCircleMetaData(M, props, values);
  else
    error('Not implemented');
  endif

endfunction

function M = Int_LoadCircleMetaData(M, props, values)
  for i=1:numel(props)
    if strcmpi(props{i}, 'type')
    elseif strcmpi(props{i}, 'uSize')
      M.uSize = values{i};
    elseif strcmpi(props{i}, 'vSize')
      M.vSize = values{i};
    elseif strcmpi(props{i}, 'uGrid')
      M.uGrid = values{i};
    elseif strcmpi(props{i}, 'vGrid')
      M.vGrid = values{i};
    elseif strcmpi(props{i}, 'uNode')
      M.uNode = values{i};
    elseif strcmpi(props{i}, 'uNodePos')
      M.uNodePos = values{i};
    else
      M = setfield(M, props{i}, values{i});
    endif
  endfor
endfunction

