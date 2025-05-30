%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022: Montpellier University / CoActions-AltRD-Emmanuel Ruffio
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
function M = HT_Model_Radiation_Init(name)
  M = struct('name', name, ...
             '__type__', 'radiationModel',  ...
             'VF', sparse(0,0),             ... % View factor matrix
             'nodes', [],                   ... % Node names associated to matrix VF
             'emissivity', [],              ... % Vector of emissivity (or matrix if they are multiple problems)
             'area', []                     ... % Vector of surface area (or matrix if they are multiple problems)
             );
endfunction
