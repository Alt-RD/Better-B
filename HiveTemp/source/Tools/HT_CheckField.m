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
%  Copyright (c) 2025 AltRD-Emmanuel Ruffio
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
% Vérification des paramètres d'entrée
function P = HT_CheckField(P, field, value, chkFunc)
  assert(nargin >= 3);
  
  if nargin < 4, chkFunc = {}; endif;
  
  % Convert to cell array
  if ~iscell(chkFunc), chkFunc = {chkFunc}; endif;
  
  lChkExist = false;
    
  for i=1:numel(chkFunc)
    if strcmpi(chkFunc{i}, 'exist')
      lChkExist = true;
      break;
    endif
  endfor
  
  if ~isfield(P, field)
    assert(~lChkExist, sprintf('Field <%s> does not exist', field));
    
    P = setfield(P, field, value);
  else  % 12/01/2023: put the function test in the "else" condition to allow invalid data to be specified as default
    for i=1:numel(chkFunc)
      if strcmp(class(chkFunc{i}), 'function_handle')
        assert(chkFunc{i}(getfield(P, field)), sprintf('Check failed: field <%s>', field));
      endif
    endfor
  endif
  
endfunction