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
% Defines the thermal emissivity of material <mat>
%
% Input arguments:
% -> mat: [object of type material] a material
% -> wavelength: [dim Mx1] [um] a matrix containing the requested wavelength
%
% Properties
% -> .model = ['constant']/'power'/'cosine' : defines the type of model used to
%                                           model the variation of emissivity with incoming light angle
% -> angle = [dim Nx1] : (default is 0) defines the zenith angle of the incoming light (in radian, 0 is vertical)
%
% Output parameters:
% v = [dim MxN] specifying the emissivity for various wavelength and various angle
% ========================================================================
%
function mat = HT_Material_SetEmissivity(mat, varargin)
  assert(HT_CheckType(mat, "material"), 'Material_SetEmissivity::Invalid material object');

  lParams = struct();

  if (numel(varargin) == 1) && isstruct(varargin{1})
    lParams = varargin{1};
    props = varargin(2:2:end);
    values = varargin(3:2:end);
  else
    props = varargin(1:2:end);
    values = varargin(2:2:end);
  endif

  assert(numel(props) == numel(values), 'Material_SetEmissivity::Invalid parameters');

  for i=1:numel(props)
    if strcmpi(props{i}, 'value')
      lParams = setfield(lParams, 'eps', values{i});
    elseif strcmpi(props{i}, 'model')
      lParams = setfield(lParams, 'epsModel', values{i});
    elseif strcmpi(props{i}, 'parameters')
      lParams = setfield(lParams, 'epsParameters', values{i});
    else
      error(sprintf('Material_SetEmissivity::Invalid argument <%s>', props{i}));
    endif
  endfor

  lParams = HT_CheckField(lParams, 'eps',             [], {@(v) (isscalar(v) && (v >= 0) && (v <= 1)) || ((columns(v) == 2) && all((v(:,2) >= 0) & (v(:,2) <= 1))  ) });
  lParams = HT_CheckField(lParams, 'epsModel',        [], {@(v) any(strcmpi(v, {'constant', 'power', 'cosine'})) });
  lParams = HT_CheckField(lParams, 'epsParameters',   [], {@(v) isrow(v) || iscolumn(v) });

  if ~isempty(lParams.eps)
    mat.eps = lParams.eps;
  endif

  if ~isempty(lParams.epsModel)
    mat.epsModel = lParams.epsModel;
  endif

  if ~isempty(lParams.epsParameters)
    mat.epsParameters = lParams.epsParameters;
  endif
endfunction

