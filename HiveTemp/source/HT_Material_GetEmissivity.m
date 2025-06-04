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
% Returns the thermal emissivity of material <mat> at specified wavelength
% <wavelength>
%
% Input arguments:
% -> mat: [object of type material] a material
% -> wavelength: [dim Mx1] [um] a matrix containing the requested wavelength
%
% Properties
% -> method = 'cubic' : defines the interpolation method used to estimate the emissivity
%                       for the specified wavelength. If no table is specified in the
%                       material object, this parameter has no effect.
% -> angle = [dim Nx1] : (default is 0) defines the zenith angle of the incoming light (in radian, 0 is vertical)
%
% Output parameters:
% v = [dim MxN] specifying the emissivity for various wavelength and various angle
% ========================================================================
%
function v = HT_Material_GetEmissivity(mat, wavelength, varargin)
  assert(HT_CheckType(mat, "material"), 'Material_GetEmissivity::Invalid material object');
  assert(all(wavelength > 0.01), 'Material_GetEmissivity::Input argument <wavelength> must be in micrometer');

  if isrow(wavelength), wavelength = wavelength'; endif

  lInterpMethod = 'cubic';
  lLightAngle = [];

  if ~isempty(varargin)
    assert(mod(numel(varargin),2) == 0, 'Material_GetEmissivity::Invalid argument count');
    prop = varargin{1:2:end};
    values = varargin{1:2:end};

    for i=1:numel(prop)
      if strcmpi(prop{i}, 'method')
        lInterpMethod = values{i};
      else
        error(sprintf('Material_GetEmissivity::Invalid argument <%s>', prop{i}));
      endif
    endfor
  endif

  if isscalar(mat.eps)
    v = repmat(mat.eps, rows(wavelength), 1);
  elseif isnumeric(mat.eps) && (columns(mat.eps) == 2)
    v = interp1(mat.eps(:,1), mat.eps(:,2), wavelength, lInterpMethod);
    if any(isna(v)) % If some wavelength lies outside the valid range, extrapolation is done using constant value
      v(wavelength < mat.eps(1,1)) = mat.eps(1,2);
      v(wavelength > mat.eps(end,1)) = mat.eps(end,2);
    endif
  elseif is_function_handle(mat.eps)
    v = mat.eps(wavelength);
  else
    error(sprintf('Material_GetEmissivity::Invalid value for emissivity <eps> in material <%s>', mat.name));
  endif

  if ~isempty(lLightAngle)
    lLightAngle((lLightAngle < 1E-10) | (lLightAngle >= pi/2)) = NA;
    lLightAngle /= pi/2;
    lLightAngle = lLightAngle(:)';

    if strcmpi(mat.epsModel, 'constant')
      v = repmat(v, 1, numel(lLightAngle));
    elseif strcmpi(mat.epsModel, 'power')
      v = v * (1 - lLightAngle.^mat.epsParameters);
    elseif strcmpi(mat.epsModel, 'cosine')
      v = v * cos(lLightAngle.^mat.epsParameters * pi/2);
    else
      error(sprintf('Material_GetEmissivity::Model <%s> is not implemented', mat.epsModel))
    endif
  endif
endfunction
