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
% Ce script affiche l'ombre projetée de la lauze sur la ruche tronc.
% Voir page 15 cahier brouillon A pour les calculs

clear variable;

pkg load matgeom;

% ==============================================================================
%                              Paramètre utilisateurs
% ==============================================================================

r = 0.4;  % [m] Rayon de la ruche tronc
h = 0.8;  % [m] Hauteur de la ruche (du cylindre)

R = 0.8;  % [m] Rayon de la lauze
H = 0.05; % [m] Hauteur de lauze

dx = 0;   % [m] Décalage de la lauze par rapport au centre du cylindre
dy = 0.3; % [m] Décalage de la lauze par rapport au centre du cylindre

alphaVec = [10 30 50]; % [0-90] Angle des rayons solaires incidents

% ==============================================================================
figure(1);
clf;
hold on;

% Draw hive body
HT_Plot_Polygon('cylinder',     'radius', 1, ...
                                'height', h/r, ...
                                'position', [0 0 -h/r], ...
                                'color', [211 171 150]/255, ...
                                'draw', true,...
                                'full', true);

% Draw hive roof
HT_Plot_Polygon('cylinder',     'radius', R/r, ...
                                'height', H/r, ...
                                'position', [0 0 0], ...
                                'color', [117 109 120]/255, ...
                                'draw', true, ...
                                'full', true);



hVec = zeros(numel(alphaVec), 1);
for i=1:numel(alphaVec)
  alpha = alphaVec(i)/360*2*pi;

  t = (-1:0.05:1)*pi/2;
  u = sin(t);
##  u = -1:0.02:1;
  y = u*r;

##  z0 = tan(alpha)*R*((r/R) - 1); % Position pour y=0;
  % Non prise en compte du décalage de la lauze par rapport au cylindre
##  z = tan(alpha)*R*(sqrt((r/R).^2-(y/R).^2) - sqrt(1-(y/R).^2));
  z = tan(alpha)*r*(sqrt(1-u.^2) - sqrt((R/r).^2-u.^2));

  % Avec décalage (dx,dy)
##  z = tan(alpha)*r*(sqrt(1-u.^2) - sqrt((R/r)^2-(u - dy/r).^2) - dx/r);

##  plot3(-sqrt(r^2-y.^2), -y, -(z-z(1)), 'displayname', sprintf('Angle=%.1f', alphaVec(i))); hold on;
  hVec(i) = plot3(sqrt(1-u.^2), u, z/r, '-+', 'linewidth', 2, 'displayname', sprintf('Angle=%.0f deg', alphaVec(i))); hold on;
endfor

hVec2 = zeros(numel(alphaVec), 1);
for i=1:numel(alphaVec)
  alpha = alphaVec(i)/360*2*pi;

  % Avec décalage (dx,dy)
  z = tan(alpha)*r*(sqrt(1-u.^2) - sqrt((R/r)^2-(u - dy/r).^2) - dx/r);
  hVec2(i) = plot3(sqrt(1-u.^2), u, z/r, 'o', 'linewidth', 4, 'displayname', sprintf('Angle(roof shift)=%.0f deg', alphaVec(i))); hold on;

  set(hVec2(i), 'color', get(hVec(i), 'color'));
endfor

xlabel('x', 'fontsize', 24, 'fontweight', 'bold');
ylabel('y', 'fontsize', 24, 'fontweight', 'bold');
zlabel('z', 'fontsize', 24, 'fontweight', 'bold');

hLegend = legend([hVec; hVec2], 'location', 'north', 'numcolumns', 2);
grid on;
axis('equal');
set(gca, 'fontsize', 16, 'fontweight', 'bold');



