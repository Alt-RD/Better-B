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
clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = './../';
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

F1 = HT_Face_Init('face1',  'axis', [1 0 0; 0 1 0]', ...
                            'norm', [0 0 1]', ...
                            'globalPosition', [0 0 0], ...
                            'size', [2;2]);
F1 = HT_Face_CreateMesh(F1, 'type', 'uniform', 'n', [5 5]);

F2 = HT_Face_Init('face1',  'axis', [1 0 0; 0 1 0]', ...
                            'norm', [0 0 1]', ...
                            'globalPosition', [-0.5 -0.5 0]', ...
                            'size', [1;1]);
F2 = HT_Face_CreateMesh(F2, 'type', 'uniform', 'n', [5 5]);


FU = HT_Face_Intersect('intersect', {F1, F2}, 'mesh', false);

figure(1);
clf;
HT_Plot_Face(F1, 'gridcolor', [0 0 1], 'displayedge', true, 'edgewidth', 2.0);
HT_Plot_Face(F2, 'gridcolor', [0 1 0], 'displayedge', true, 'edgewidth', 2.0);
HT_Plot_Face(FU, 'gridcolor', [1 0 0], 'displayedge', true, 'edgewidth', 2.0, 'edgecolor', 'red');
xlabel('x');
ylabel('y');
zlabel('z');
return;


