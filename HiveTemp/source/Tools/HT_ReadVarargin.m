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
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2023-2025: CoActions-AltRD-Emmanuel Ruffio
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
function [lKeys lValues lParams lOptions] = HT_ReadVarargin(args)
  lKeys = {};
  lValues = {};
  lParams = struct();
  lOptions = struct();

  if isempty(args)
    return;
  endif

  if isstruct(args{1})        % If first parameter is a struct, parameters are loaded from it
    lParams = args{1};

    if (mod(numel(args), 2) == 0)
      assert(isstruct(args{end}), 'Invalid arguments <args>. End param must be a struct');

      lOptions = args{end};
      lKeys = args(2:2:(end-1));
      lValues = args(3:2:(end-1));
    else
      lKeys = args(2:2:end);
      lValues = args(3:2:end);
    endif
  else
    if (mod(numel(args), 2) == 1)
      assert(isstruct(args{end}), 'Invalid arguments <args>. End param must be a struct');

      lOptions = args{end};
      lKeys = args(1:2:(end-1));
      lValues = args(2:2:(end-1));
    else
      lKeys = args(1:2:end);
      lValues = args(2:2:end);
    endif
  endif

  if isargout(3)
    lFieldNames = fieldnames(lParams);

    % Add to the <lParams> structure all parameters specified
    for i=1:numel(lKeys)
      lParams = setfield(lParams, lKeys{i}, lValues{i});
      if ~any(strcmp(lFieldNames, lKeys{i}))
        lFieldNames = [lFieldNames; lKeys{i}];
      endif
    endfor
  endif
endfunction
