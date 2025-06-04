%  This file is part of project HiveModel.
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

clear variable;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../../../HiveTemp/source/';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

% ==============================================================================
%                              Paramètre utilisateurs
% ==============================================================================

r = 0.4;  % [m] Rayon de la ruche tronc
h = 0.8;  % [m] Hauteur de la ruche (du cylindre)
nh = 10;
nTheta = 16;

R = 0.8;  % [m] Rayon de la lauze
H = 0.05; % [m] Hauteur de lauze

dx = 0;   % [m] Décalage de la lauze par rapport au centre du cylindre
dy = 0;   % [m] Décalage de la lauze par rapport au centre du cylindre

alphaVec = [10 30]; % [0-90] Angle des rayons solaires incidents
azimutVec = [ 0 0];

plotRoof = true;
plotShiftedRoof = true;

materialEmissivityType = 'power';
materialEmissivityPower = 4;

% ==============================================================================

lAngle = [azimutVec; alphaVec]/360*2*pi;
heightVec = (0:(nh))'/(nh) * h;
thetaVec = (0:nTheta)'/nTheta*2*pi;

[A J S] = HT_Shadow_GetDiskOverCylinder(R, r, heightVec, thetaVec, lAngle, ...
                      'deltaDisk', [dx dy], ...
                      'epsType', materialEmissivityType, ...
                      'power', materialEmissivityPower);
X = reshape(J(:,1), nTheta, nh)';

% ==============================================================================
%   Plot the 3D cylinder showing the ratio of shadowed area over total area
% ==============================================================================
if true
figure(1);
clf;

for i=1:min(numel(alphaVec), 4)
  subplot(2,numel(alphaVec),i);
  [V F] = HT_Plot_Polygon('cylinder',   'radius', r, ...
                                        'n', nTheta, ...
                                        'height', heightVec, ...
                                        'position', [0 0 0], ...
                                        'draw', true,...
                                        'sort', 'theta', ...
                                        'full', false, ...
                                        'color', A(:,i) );

  if plotRoof
  [V F] = HT_Plot_Polygon('cylinder',   'radius', R, ...
                                        'n', nTheta, ...
                                        'height', H, ...
                                        'position', [0 0 heightVec(end)], ...
                                        'color', [117 109 120]/255, ...
                                        'draw', true,...
                                        'full', true );
  endif

  if any([dx dy] != 0) && plotShiftedRoof
  HT_Plot_Polygon('disk',         'radius', R, ...
                                  'position', [dx dy heightVec(end)], ...
                                  'color', 'none', ...
                                  'draw', true);
  endif

  xlabel('x', 'fontsize', 24, 'fontweight', 'bold');
  ylabel('y', 'fontsize', 24, 'fontweight', 'bold');
  zlabel('z', 'fontsize', 24, 'fontweight', 'bold');

  grid on;
  axis('equal');
  set(gca, 'fontsize', 16, 'fontweight', 'bold');
  title(sprintf("Shadow fraction\nElevation=%.0f deg, Azimuth=%.0f deg", alphaVec(i), azimutVec(i)));

##  campos( [14.7963   -4.1294    5.7743]);
##  camtarget( [0 0.1 0.4250 ]);
##  camup([0 0 1]);

  Jdensity = J(:,i) ./ S;

  subplot(2,numel(alphaVec),i+numel(alphaVec));
  [V F] = HT_Plot_Polygon('cylinder',   'radius', r, ...
                                        'n', nTheta, ...
                                        'height', heightVec, ...
                                        'position', [0 0 0], ...
                                        'draw', true,...
                                        'sort', 'theta', ...
                                        'full', false, ...
                                        'color', Jdensity / max(max(Jdensity)));

  if plotRoof
    [V F] = HT_Plot_Polygon('cylinder',   'radius', R, ...
                                          'n', nTheta, ...
                                          'height', H, ...
                                          'position', [0 0 heightVec(end)], ...
                                          'color', [117 109 120]/255, ...
                                          'draw', true,...
                                          'full', true );
  endif

  if any([dx dy] != 0) && plotShiftedRoof
    HT_Plot_Polygon('disk',         'radius', R, ...
                                    'position', [dx dy heightVec(end)], ...
                                    'color', 'none', ...
                                    'draw', true);
  endif

  xlabel('x', 'fontsize', 24, 'fontweight', 'bold');
  ylabel('y', 'fontsize', 24, 'fontweight', 'bold');
  zlabel('z', 'fontsize', 24, 'fontweight', 'bold');

  grid on;
  axis('equal');
  set(gca, 'fontsize', 16, 'fontweight', 'bold');
  title(sprintf("Irrradiance\nElevation=%.0f deg, Azimuth=%.0f deg", alphaVec(i), azimutVec(i)));
##
##  campos( [14.7963   -4.1294    5.7743]);
##  camtarget( [0 0.1 0.4250 ]);
##  camup([0 0 1]);
endfor
endif


% ==============================================================================
%   Plot the total energy absorbed by the hive with respect to the angle
% ==============================================================================
figure(2);
clf;

lParams = struct('epsType', {'constant', 'power', 'power', 'power', 'cosine', 'cosine'}, ...
                 'power',   {0         , 2      , 4      , 8     , 2       , 4}, ...
                 'displayname', '%s.%d');

lAlphaVec = (0:0.02:1).^0.8*90;
lAzimuth = zeros(size(lAlphaVec));

for i=1:numel(lParams)
  [A J S] = HT_Shadow_GetDiskOverCylinder(R, r, [heightVec(1) heightVec(end)], [], ...
                        [lAzimuth; lAlphaVec]/360*2*pi, ...
                        'deltaDisk', [0 0], ...
                        'epsType', lParams(i).epsType, ...
                        'power', lParams(i).power);

  plot(lAlphaVec, J,  'displayname', sprintf(lParams(i).displayname, lParams(i).epsType, lParams(i).power), ...
                      'linewidth', 2); hold on;
endfor

legend();

grid on;
set(gca, 'fontsize', 16, 'fontweight', 'bold');
xlabel('Elevation angle [deg]', 'fontsize', 20, 'fontweight', 'bold');
ylabel('Absorbed energy over solar irradiance [m^2]', 'fontsize', 20, 'fontweight', 'bold');
title("Energy absorbed by direct illumination", 'fontsize', 20, 'fontweight', 'bold');

return;

