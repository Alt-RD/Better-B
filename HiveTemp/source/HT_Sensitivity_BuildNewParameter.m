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
%  Copyright (c) 2022-2025 AltRD: Emmanuel Ruffio
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
## Three different model for angle dependent emissivity
##  The <angle> is the normalized angle between the incoming light and the
##  surface normal. <angle> (in [0;1]) is then the angle divided by pi/2.
##  'constant':
##  'power': eps(angle) = eps * (1 - angle.^mat.epsParameters);
##  'cosine': eps(angle) = eps * cos(angle.^mat.epsParameters * pi/2);
##
function lParameters = HT_Sensitivity_BuildNewParameter(varargin)
  assert(numel(varargin) > 0);

  lParameters = struct('name', '', ...
                       'displayName', '', ...
                       'fields', [], ...
                       'epsilon', 1E-3, ...
                       'type', 'scale', ...
                       'setterNext', @(object, params) setfield(params, object.fields{:}, object.setter(getfield(params, object.fields{:}), object.getter(getfield(params, object.fields{:}))*(1+object.epsilon))), ...
                       'setterPrev', @(object, params) setfield(params, object.fields{:}, object.setter(getfield(params, object.fields{:}), object.getter(getfield(params, object.fields{:}))*(1-object.epsilon))), ...
                       'getter', @(fieldObject) fieldObject, ... % Convert field object to value. For example: material lambda -> lambda value. By default, field are assumed to be standard scalar value
                       'setter', @(fieldObject, value) value, ... % By default, field are assumed to be standard scalar value
                       'compute', @(object, params, Xp, Xn) (Xn-Xp) / (2*object.getter(getfield(params, object.fields{:}))*object.epsilon));

  lOptions = struct('verbose', true);

  [props, values, lParameters, lOptions] = HT_ReadVarargin(lParameters, lOptions, varargin, 'struct array');

  for i=1:numel(lParameters)
    lParam = lParameters(i);
    HT_CheckField(lParam, 'name',              '', {'exist', @(v) ischar(v)});
    HT_CheckField(lParam, 'displayName',       '', {@(v) ischar(v)});
    HT_CheckField(lParam, 'fields',            '', {@(v) ischar(v) || iscellstr(v) });
    HT_CheckField(lParam, 'epsilon',           1E-3, {@(v) isnumeric(v) && isscalar(v) });
    HT_CheckField(lParam, 'type',              'scale', {@(v) ischar(v) && any(strcmpi(v, {'scale', 'add'})) });

    % Convert fields string to cellstr
    lParameters(i).fields = strsplit(lParameters(i).fields, '.');

##    lGetterInfo = functions(lParameters(i).getter);
##    lParamList = lGetterInfo.function(strchr(lGetterInfo.function, '(', 1):strchr(lGetterInfo.function, ')', 1));
##    lParamList = strsplit(lParamList, ',');
##    if numel(lParamList) == 1
##      lParameters(i).getter = @(object, params) lParameters(i).getter(getfield(params, object.fields{:}));
##    endif
  endfor


endfunction


