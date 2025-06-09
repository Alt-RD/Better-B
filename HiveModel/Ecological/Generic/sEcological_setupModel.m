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

% ========================================================================
%                            HIVE THERMAL MODEL
% ========================================================================

% Create the figure(1) that will contain the geometry
if lPlotOptions.plotGeometry
  if lPlotOptions.onlyOneFigure, figure(1); clf; subplot(1,1,1); ...
  else figure(1); clf; endif
  set(gca, 'FontSize', 14);
  xlabel('X', 'fontsize', 16, 'fontweight', 'bold');
  ylabel('Y', 'fontsize', 16, 'fontweight', 'bold');
  zlabel('Z', 'fontsize', 16, 'fontweight', 'bold');
  axis equal;
endif

% Create model main object with name "main"
lModel = HT_Model_Init('', 'main');

% Create the node corresponding to the outside air
% It will be connected to face when necessary along the model building process
lAirExtNode = HT_Node_Init('nAirExt', 'mode', 'distributed');

% ========================================================================
% Build air body part
% ========================================================================
% Yet, HiveTemp library does not allow to define 1D conduction in custom shape
% Only rectangle shape is possible with Conduction1D. It could be extended in
% the future.

% Internal radius of the hive computed between two parallel faces
lIntRadius = lHiveParams.extRadius - sum([lHiveParams.wall.intThickness + lHiveParams.wall.insulationThickness + lHiveParams.wall.extThickness]);
lIntEdgeLength = 2*tan(pi/8)*lIntRadius;

lIntArea = 8*lIntRadius^2*(sqrt(2)-1);  % Area of the internal plane of the body
lIntSimplifiedArea = lIntEdgeLength*lIntEdgeLength;
lAirCapacityCorrectionFactor = lIntArea / lIntSimplifiedArea;
lAirCorrected = HT_Material_SetCapacity(lHiveMat.internalAir, lAirCapacityCorrectionFactor*HT_Material_GetCapacity(lHiveMat.internalAir));

lAverageRoofHeight = tan(lHiveParams.roofTilt) * lIntRadius/2;

[lMod_AirCylinder, lMod_AirCylinderFaces] = HT_Model_Conduction1D("mAirBody", ...
                                          struct( "length", {{ lHiveParams.wall.baseThickness; ...
                                                              lHiveParams.bodyHeight; ...
                                                              lHiveParams.overframeThickness; ...
                                                              lHiveParams.ovfHeight; ...
                                                              lAverageRoofHeight; ...
                                                              lHiveParams.roofThickness}}, ... Hive height
                                                  "material", [ lHiveMat.walls; ...
                                                                lAirCorrected;...
                                                                lHiveMat.overframe; ...
                                                                lAirCorrected; ...
                                                                lAirCorrected; ...
                                                                lHiveMat.roof], ...
                                                  "rc", [ 1/lConvection.body_bot; ...
                                                          1/lConvection.body_top; ...
                                                          1/lConvection.ovf_bot; ...
                                                          1/lConvection.ovf_top; ...
                                                          0], ...
                                                  "dim", lIntEdgeLength*ones(2,1), ...
                                                  "axis", [0 1 0;
                                                           0 0 1;
                                                           1 0 0], ...
                                                  "globalPosition", [-lIntEdgeLength;-lIntEdgeLength;0]/2, ...
                                                  "n", {{  lMesh.wall.bottom; ...
                                                          lMesh.bodyHeight; ...
                                                          lMesh.wall.ovf; ...
                                                          lMesh.ovfHeight; ...
                                                          1;...
                                                          lMesh.wall.roof}}, ...
                                                  "mergeFaces", false, ...
                                                  "gridType", {{'hf'; 'hh'; 'hh'; 'hf'; 'ff'; 'hh' }}), lOptions);

clear lIntArea lAverageRoofHeight;

if lPlotOptions.plotGeometry
  HT_Plot_Face(lMod_AirCylinderFaces, l3DPlotOptions); % Display the geometry
endif

lModel = HT_Model_Merge(lModel, lMod_AirCylinder, lOptions);

lFaceList = HT_Face_Find(lMod_AirCylinderFaces, 'orientation', [1 0 0]);
lFaceUnion = HT_Face_Union("fSideXp", lFaceList, 'mesh', true, 'coplanar', true, 'allowIntersect', false);
##lFaceUnion = HT_Face_Transform(lFaceUnion, "resize", ...
##              struct('size', lIntEdgeLength, ...
##                      'direction', [0; 1; 0], ...
##                      'center', [0.5 0.5]));
lFaceUnion = HT_Face_Transform(lFaceUnion, "position", ...
              struct('position', [lIntRadius; NA; NA], ...
                     'axis', [0, 0;
                              1, 0;
                              0, 1], ...
                     'center', [0.5; 0]));

% ========================================================================
% Build the hive sides (8 sides)
% ========================================================================
% Define the face that will be used to define thermal connection
lIntBaseFace =        lMod_AirCylinderFaces(4);
lIntAirFace =         lMod_AirCylinderFaces(10);
lIntOverframeFace =   lMod_AirCylinderFaces(16);
lIntAirOverframeFace = lMod_AirCylinderFaces(22);
lIntAirRoofFace =     lMod_AirCylinderFaces(28);
lIntRoofFace =        lMod_AirCylinderFaces(34);

% Store the external faces for later use (in setupCommand file)
lLateralOutsideFaces = cell(8,1);
lRoofOutsideFaces = cell(8,1);

for i=1:8
  lFaceBase = HT_Face_Transform(lFaceUnion, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));

  [lModSide, lModSide_faces] = HT_Model_Conduction3D(...
    sprintf('mSide_%d', i),  ...
    struct( 'base',           lFaceBase,                                ... % Internal face
            'length',         [lHiveParams.wall.intThickness; ...
                               lHiveParams.wall.insulationThickness; ...
                               lHiveParams.wall.extThickness],   ...
            'material',       [lHiveMat.walls; ...
                               lHiveMat.wallInsulation; ...
                               lHiveMat.walls], ...                    ... % Retrieve the material specified in the header (generally air)
            'n',              [NA NA lMesh.wall.intLayer ...
                                     lMesh.wall.insulationLayer ...
                                     lMesh.wall.extLayer],                    ...
            'zGridType',      'hh'                                     ... % Specify full node size on boundaries (instead of half). This mode must be selected to allow 1 node only in the z-direction
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lModSide, lOptions);

  lLateralOutsideFaces{i} = lModSide_faces(FaceZP);
  lRoofOutsideFaces{i} = lModSide_faces(FaceYP);
  % List of faces <lModSide_faces> are in the following order:
  % Xm, Xp, Ym, Yp, Zm, Zp

  % Connection to the external air (starting from bottom to top)
  % Connect this side to the external air (bottom face)
  assert(norm(lModSide_faces(FaceYM).norm + lHiveParams.globalAxis(:,Z), 1) < 1E-12);

  lMod_Tmp = HT_Model_Connect(sprintf('mcSide_%d_BotAir', i), ...
                                  lAirExtNode,                      ... % Node of outside air
                                  lModSide_faces(FaceZM),           ... % Top face of the roof
                                  'g', lConvection.ext_bot,       ... % Convection coefficient specified by user
                                  lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect the laterial side to the external air
  lMod_Tmp = HT_Model_Connect(sprintf('mcSide_%d_SideAir', i), ...
               lAirExtNode,                     ... % Node of outside air
               lModSide_faces(FaceZM),          ... % Top face of the roof
               'g', lConvection.ext_walls,      ... % Convection coefficient specified by user
               lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the external air (top face)
  assert(norm(lModSide_faces(FaceYP).norm - lHiveParams.globalAxis(:,Z), 1) < 1E-12);

  lMod_Tmp = HT_Model_Connect(sprintf('mcSide_%d_TopAir', i), ...
                                  lAirExtNode,                      ... % Node of outside air
                                  lModSide_faces(FaceZP),           ... % Top face of the roof
                                  'g', lConvection.roof_top,       ... % Convection coefficient specified by user
                                  lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connection to the internal parts (starting from bottom)
  % Connect this side to the hive base
  lIntFaceRotated = HT_Face_Transform(lIntBaseFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_Base', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the hive air body
  lIntFaceRotated = HT_Face_Transform(lIntAirFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_BodyAir', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the overframe wood layer
  lIntFaceRotated = HT_Face_Transform(lIntOverframeFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_Ovf', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the hive air overframe
  lIntFaceRotated = HT_Face_Transform(lIntAirOverframeFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_OvfAir', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the hive air roof
  lIntFaceRotated = HT_Face_Transform(lIntAirRoofFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_RoofAir', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect this side to the hive roof
  lIntFaceRotated = HT_Face_Transform(lIntRoofFace, "rotate", struct('angle', -(i-1)*2*pi/8, ...
                                                               'axis', 'Z', ...
                                                               'center', [0;0;0]));
  lMod_Tmp = HT_Model_ConnectFaces(sprintf('mcSide_%d_Roof', i),         ...
                                  lIntFaceRotated,            ... % External face of the internal body cylinder
                                  lModSide_faces(FaceZM),         ... % Internal face of the wall
                                  struct('intersect', true),                       ... % Convection coefficient specified by user
                                  setfield(lOptions, 'checkGeometry', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  if lPlotOptions.plotGeometry
    HT_Plot_Face(lModSide_faces, l3DPlotOptions); % Display the geometry
  endif

endfor

disp('==============================================');
disp(sprintf('Model built in %.2f s', lModel.timer.buildingTime));


