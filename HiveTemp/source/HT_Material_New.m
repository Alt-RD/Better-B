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
function M = HT_Material_New(varargin)
  assert(numel(varargin) > 0);
  assert(ischar(varargin{1}), 'Invalid argument. Should be string');

  lParamList = {'name', 'rho', 'cp', 'eps', 'lambda', 'color'};
  lLoadFromName = ~any(strcmpi(varargin{1}, lParamList));

  if lLoadFromName   % Material name specified
    M = INT_GetFromName(varargin{1});
    prop = varargin(2:2:end);
    value = varargin(3:2:end);
  else
    M = struct('__type__', 'material');
    prop = varargin(1:2:end);
    value = varargin(2:2:end);
  endif

  assert(numel(prop) == numel(value), 'Invalid input parameters'); % Even number of arguments

  for i = 1:numel(prop)
    lParamName = prop{i};
    assert(any(strcmpi(prop{i}, lParamList)), sprintf('Invalid parameter <%s>', prop{i}));
    M = setfield(M, prop{i}, value{i});
  endfor

  M = HT_CheckField(M, 'color',             [0, 0, 0]);
  M = HT_CheckField(M, 'eps',               1.0);
  M = HT_CheckField(M, 'epsModel',          'constant', {@(v) any(strcmpi(v, {'constant', 'power', 'cosine'}))});
  M = HT_CheckField(M, 'epsParameters',     1.0); % This field holds the parameters used by the epsModel
  % Mandatory fields
  M = HT_CheckField(M, 'rho',   0,          {'exist'});
  M = HT_CheckField(M, 'cp',    0,          {'exist'});
  M = HT_CheckField(M, 'lambda',   0,       {'exist'});

  % For hygro: parameters with dimension [1x3] handle linear temperature and water density dependency
  M = HT_CheckField(M, 'molecularDiffusion',   NA, {@(v) isna(v) || all(size(v) == [1 3])});
  M = HT_CheckField(M, 'liquidPermeability',   NA, {@(v) isna(v) || all(size(v) == [1 3])});
  M = HT_CheckField(M, 'molarMass',            NA, {@(v) isscalar(v) && isnumeric(v)});
  M = HT_CheckField(M, 'liquidRho',            NA, {@(v) isscalar(v) && isnumeric(v)});
  M = HT_CheckField(M, 'Tref',                 293.15, {@(v) isscalar(v) && isnumeric(v)}); % Reference temperature used for temperature dependant parameter
  M = HT_CheckField(M, 'sorptionCurve',        NA, {@(v) isna(v) || is_function_handle(v) || isstruct(v)}); % Reference temperature used for temperature dependant parameter
  % sorptionCurve can be a function handle f(rho/rhoanh,T) or f(phi,T)
  % Or
  %             -> cell array of function_handle or function_handle if only one layer or similar layer curve
  %             -> function_handle: [phi DphiDrho DphiDT] = @(rhoAdim, T)
  %                 if DphiDrho or DphiDT, they are computed by finite difference
  %             -> 2d interpolation array with structure of size (1) of size = nLayer
  %                   "data" : [LxC] matrix containing relative humidity values
  %                   "rho"  : [Lx1] scaled rho values (by rho_anhydre) of each line of the matrix
  %                   "T"    : [1xC] temperature values of each column of the matrix

endfunction

function M = INT_GetFromName(name)
  HT_ImportConstants();

  if      strcmpi(name, 'wood')
    M = HT_Material_New(  'name', name, ...
                          'color', [106 77 21]/255, ...
                          'rho', 600, ...
                          'cp', 1600, ...
                          'lambda', 0.16, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5; ...
                                  HT_WAVELENGTH_IR_MED   0.9]);
  elseif  strcmpi(name, 'polystyrene')
    M = HT_Material_New(  'name', name, ...
                          'color', [192 192 192]/255, ...
                          'rho', 20, ...  30 à 40kg/m3
                          'cp', 1450, ...
                          'lambda', 0.035, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.9]);
  elseif  strcmpi(name, 'air')
    M = HT_Material_New(  'name', name, ...
                          'color', [3  176 226]/255, ...
                          'rho', 1.2, ...
                          'cp', 1009, ...
                          'lambda', 0.025, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.9]);
  elseif  strcmpi(name, 'steal')
    M = HT_Material_New(  'name', name, ...
                          'color', [117 109 120]/255, ...
                          'rho', 7800, ...
                          'cp', 450 , ...
                          'lambda', 50, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.45;
                                  HT_WAVELENGTH_IR_MED   0.9]);
  elseif  strcmpi(name, 'dry ground')
    M = HT_Material_New(  'name', name, ...
                          'color', [185 122 87]/255, ...
                          'rho', 1200, ...
                          'cp', 1250 , ...
                          'lambda', 0.7, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5;
                                  HT_WAVELENGTH_IR_MED   0.95]);
  elseif  strcmpi(name, 'wet ground')
    M = HT_Material_New(  'name', name, ...
                          'color', [185 122 87]/255, ...
                          'rho', 1700, ...
                          'cp', 1250 , ...
                          'lambda', 1.7, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5;
                                  HT_WAVELENGTH_IR_MED   0.95]);
  elseif  strcmpi(name, 'granit') % Wikipedia
    M = HT_Material_New(  'name', name, ...
                          'color', [125 81 81]/255, ...
                          'rho', 2500, ...
                          'cp', 790 , ...
                          'lambda', 2.2, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5;
                                  HT_WAVELENGTH_IR_MED   0.95]);
  elseif  strcmpi(name, 'brick') % https://material-properties.org/fr/brique-densite-capacite-thermique-conductivite-thermique/?utm_content=cmp-true
    M = HT_Material_New(  'name', name, ...
                          'color', [255 129 39]/255, ...
                          'rho', 1700, ...
                          'cp', 800 , ...
                          'lambda', 1.3, ...
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5;
                                  HT_WAVELENGTH_IR_MED   0.95]);
  elseif  strcmpi(name, 'concrete') % https://material-properties.org/fr/brique-densite-capacite-thermique-conductivite-thermique/?utm_content=cmp-true
    M = HT_Material_New(  'name', name, ...
                          'color', [195 195 195]/255, ...
                          'rho', 2200, ...
                          'cp', 880 , ...     % Wikipedia
                          'lambda', 0.92, ... % Wikipedia
                          'eps', [HT_WAVELENGTH_VISIBLE  0.5;
                                  HT_WAVELENGTH_IR_MED   0.95]);
  else
    error(sprintf('Unknown material <%s>', name));
  endif
endfunction

