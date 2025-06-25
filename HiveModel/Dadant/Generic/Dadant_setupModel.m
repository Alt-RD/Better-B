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
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2023-2025: CoActions-AltRD-Emmanuel Ruffio
%  Author: emmanuel.ruffio@gmail.com
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
##function M = HiveModel_Dadant_setupCommands(hiveOptions, hiveMat, convection, contact, radiation, mesh, plotOptions, plot3dOptions, plotColors, options)
function [M faces volumes] = Dadant_setupModel(varargin)
  HT_ImportConstants();
  assert(nargin > 0);

  if numel(varargin) == 1 && isstruct(varargin{1})
    lParams = varargin{1};

    lHiveParams =   lParams.hiveParams;
    lHiveMat =      lParams.hiveMat;
    lConvection =   lParams.convection;
    lContact =      lParams.contact;
    lRadiation =    lParams.radiation;
    lMesh =         lParams.mesh;
    lOptions =      lParams.options;
    lPlotOptions  = lParams.plotOptions;
    lPlotColors =   lParams.plotColors;
    l3DPlotOptions = lParams.plotOptions3d;
  else
    error('Invalid parameter');  prop = varargin(1:2:end);
 ##    prop = varargin(1:2:end);
##    value = varargin(2:2:end);
##    assert(numel(prop) == numel(value), 'Invalid properties');
##
##    lHiveParams = struct();
##    lHiveMat = struct();
##    lConvection = struct();
##    lContact = struct();
##    lRadiation = struct();
##    lMesh = struct();
##    lPlotOptions = struct();
##    l3DPlotOptions = struct();
##    lPlotColors = struct();
##    lOptions = struct();
##
##    for i=1:numel(prop)
##      if strcmpi(prop{i}, 'hiveParams'), lHiveParams = value{i};
##      elseif strcmpi(prop{i}, 'hiveMaterial'), lHiveMat = value{i};
##      elseif strcmpi(prop{i}, 'convection'), lConvection = value{i};
##      elseif strcmpi(prop{i}, 'contact'), lContact = value{i};
##      elseif strcmpi(prop{i}, 'radiation'), lRadiation = value{i};
##      elseif strcmpi(prop{i}, 'mesh'), lMesh = value{i};
##      elseif strcmpi(prop{i}, 'plotOptions'), lPlotOptions = value{i};
##      elseif strcmpi(prop{i}, '3dPlotOptions'), l3DPlotOptions = value{i};
##      elseif strcmpi(prop{i}, 'plotColors'), lPlotColors = value{i};
##      elseif strcmpi(prop{i}, 'options'), lOptions = value{i};
##      else
##        error(sprintf('Invalid property <%s>', prop{i}));
##      endif
##    endfor
  endif

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


  % ============================== Hivebody 1 2 ==================================
  % Create the hive body with dimension 0.5x0.5x0.5m
  [lMod_HiveBody, lMod_HiveBody_faces, ~, lMod_HiveBodyVolumes] = HT_Model_Conduction1D(            ... % Build a two layers 1D heat conduction problem
                            "mHB",                                        ... % Model name of the hive body part 1 and 2
                            'material', lHiveMat.bodyAir,                 ... % Air material for the hive body parts
                            'dim', lHiveParams.dimension,                 ... % Dimensions along X and Y direction
                            'axis', [0 0 1; 1 0 0; 0 1 0]',               ... % Local axis (x,y,z)=(Z,X,Y) so that x is oriented towards Z+
                            'length', {lHiveParams.bodyHeight(1); lHiveParams.bodyHeight(2)},... % Length of the two layers (part1 and part2)
                            'n', {lMesh.hivebody(1); lMesh.hivebody(2)},  ... % Number of nodes in each layer
                            'mergeFaces', false,...
                            'gridType', {'hf'; 'fh'}, ...
                            'boundaryNames', {{ 'nHBBot', '' }; {'', 'nHBTop'}}, ...
                            lOptions);
  % Faces are returned following local axis. Faces are permuted to match global axis
  lMod_HiveBody1_faces = lMod_HiveBody_faces(    [FaceYM; FaceYP; FaceZM; FaceZP; FaceXM; FaceXP]);
  lMod_HiveBody2_faces = lMod_HiveBody_faces(6 + [FaceYM; FaceYP; FaceZM; FaceZP; FaceXM; FaceXP]);
  clear lMod_HiveBody_faces;

  % Merge the 0D Model to the main model
  lModel = HT_Model_Merge(lModel, lMod_HiveBody, lOptions);

  if lPlotOptions.plotGeometry
  HT_Plot_Face(lMod_HiveBody1_faces, l3DPlotOptions); % Display the geometry
  ##HT_Plot_Face(lMod_HiveBody2_faces, l3DPlotOptions); % Display the geometry
  endif

  % ================================= Sidewalls ==================================
  % Build hive side walls 1 and 2(XM,XP,YM,YP)
  lMod_SideWall1_facesCell = cell(4,1);   % Backup the faces of sidewalls for later use
  lMod_SideWall2_facesCell = cell(4,1);   % Backup the faces of sidewalls for later use

  for i=1:4
    [lMod_SideWall1, lMod_SideWall1_faces] = HT_Model_Conduction3D(...
      sprintf('mSW1_%s', FaceDIRstrL{i}), ...
      struct( 'base',     lMod_HiveBody1_faces(i),    ...
              'length',   lHiveParams.wallThickness,  ...
              'material', lHiveMat.sidewalls,         ...
              'n',        [NA NA lMesh.sidewalls]     ... % NA is used since mesh is extracted from face <base>
              ),...
      lOptions);

    lModel = HT_Model_Merge(lModel, lMod_SideWall1, lOptions);
    lMod_SideWall1_facesCell{i} = lMod_SideWall1_faces;

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_SideWall1_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Connect sidewalls 1 to the hive body 1
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcBody1_SideWall1_%s', FaceDIRstrL{i}), ...
      lMod_HiveBody1_faces(i), ...
      lMod_SideWall1_faces(FaceZM),...
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect sidewalls 1 to external air
    lMod_Tmp = HT_Model_Connect(sprintf('mcSideWall1_%s_ExtAir', FaceDIRstrL{i}), ...
                 lAirExtNode,                     ... % Node of outside air
                 lMod_SideWall1_faces(FaceZP),         ... % Top face of the roof
                 'g', lConvection.walls1_ext{i},  ... % Convection coefficient specified by user
                 lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Hive Body part 2
    [lMod_SideWall2, lMod_SideWall2_faces] = HT_Model_Conduction3D(...
      sprintf('mSW2_%s', FaceDIRstrL{i}), ...
      struct( 'base',     lMod_HiveBody2_faces(i),    ...
              'length',   lHiveParams.wallThickness,  ...
              'material', lHiveMat.sidewalls,         ...
              'n',        [NA NA lMesh.sidewalls]     ... % NA is used since mesh is extracted from face <base>
              ),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_SideWall2, lOptions);
    lMod_SideWall2_facesCell{i} = lMod_SideWall2_faces;

    % Connect sidewalls 2 to the hive body 2
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcBody2_SideWal2_%s', FaceDIRstrL{i}), ...
      lMod_HiveBody2_faces(i), ...
      lMod_SideWall2_faces(FaceZM),...
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_SideWall2_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Find the face among "hiveBody1 faces" that is oriented to Z+
    lFaceZp = HT_Face_Find(lMod_SideWall1_faces, 'orientation', lHiveParams.globalAxis(:,Z), 'numberCheck', 1, 'returnType', true);
    % Find the face among "hiveBody2 faces" that is oriented to Z-
    lFaceZm = HT_Face_Find(lMod_SideWall2_faces, 'orientation', -lHiveParams.globalAxis(:,Z), 'numberCheck', 1, 'returnType', true);

    % Create a thermal connection between them
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcHiveBody12_%s', FaceDIRstrL{i}),   ... % Connection model name
      lFaceZp,...
      lFaceZm,...
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  endfor

  [lMod_SideWallBot, lMod_SideWallBot_faces] = HT_Model_Conduction3D(...
    'mSWlBot', ...
    struct( 'base',     lMod_HiveBody1_faces(FaceZM),    ...
            'length',   lHiveParams.wallThickness,  ...
            'material', lHiveMat.sidewalls,         ...
            'n',        [NA NA lMesh.sidewalls],    ... % NA is used since mesh is extracted from face <base>
            'zGridType', 'ff' ...
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_SideWallBot, lOptions);

  % Reorder the face to match global axis
  lMod_SideWallBot_faces = HT_Face_Find(lMod_SideWallBot_faces, 'orientation', FaceDIRmat, 'numberCheck', 1, 'returnType', true);


  % Plot geometry
  if lPlotOptions.plotGeometry
    HT_Plot_Face(lMod_SideWallBot_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
  endif

  % Connect sidewalls 2 to the hive body 2
  lMod_Tmp = HT_Model_ConnectFaces(...
    'mcHiveBody1WallBot', ...
    lMod_SideWallBot_faces{FaceZP}, ...
    lMod_HiveBody1_faces(FaceZM),...
    struct(),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % ================================= Overframe ==================================
  [lMod_Overframe, lMod_Overframe_faces] = HT_Model_Conduction3D(...
    'mOvf', ...
    struct( 'base',     lMod_HiveBody2_faces(FaceZP),     ...
            'length',   lHiveParams.overframeThickness,   ...
            'material', lHiveMat.overframe,               ...
            'n',        lMesh.overframe,                  ...
            'zGridType','ff',                             ...
            'xGridType','ff',                             ...
            'yGridType','ff'...
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Overframe, lOptions);

  % Reorder the face to match global axis
  % lMod_Overframe_faces is now a cell array instead of a struct array
  lMod_Overframe_faces = HT_Face_Find(lMod_Overframe_faces, 'orientation', FaceDIRmat, 'numberCheck', 1, 'returnType', true);

  % Plot geometry
  if lPlotOptions.plotGeometry
    HT_Plot_Face(lMod_Overframe_faces, l3DPlotOptions);
  endif

  % Connect overframe to the hive body 2
  lMod_Tmp = HT_Model_ConnectFaces(...
    'mcOvframeHiveBody2', ...
    lMod_HiveBody2_faces(FaceZP), ...
    lMod_Overframe_faces{FaceZM},...
    struct(),...
    setfield(lOptions, 'strictMatch', false));
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  lMod_OvfSideWall_facesCell = cell(4,1);    % Backup faces for later use
  % Add surrounding wood sides
  for i=1:4
    % Retrieve the base face of overframe layer

    [lMod_OvfSW, lMod_OvfSW_faces] = HT_Model_Conduction3D(...
      sprintf('mOvfSW_%s', FaceDIRstrL{i}),           ... % Model name
      struct( 'base',     lMod_Overframe_faces{i},    ... % Base face over which the 1D model is built
              'length',   lHiveParams.wallThickness,  ...
              'material', lHiveMat.sidewalls,         ...
              'n',        [NA NA lMesh.sidewalls],    ... % NA is used since mesh is extracted from face <base>
              'zGridType', 'hh',...                   ... % Half nodes are used on boundaries
              'xGridType', 'ff',...
              'yGridType', 'ff'...
              ),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_OvfSW, lOptions);
    lMod_OvfSideWall_facesCell{i} = lMod_OvfSW_faces;

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_OvfSW_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Connect overframe sidewalls to the overframe
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcOvfSW_%s', FaceDIRstrL{i}),  ...
      lMod_Overframe_faces{i},                ...
      lMod_OvfSW_faces(FaceZM),               ... % For 1D conduction problem, FaceXM is always the base face
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect overframe sidewalls to hivebody 2 sidewalls
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcOvfSW_HB2SW_%s', FaceDIRstrL{i}),  ...
      HT_Face_Find(lMod_SideWall2_facesCell{i}, 'orientation', lHiveParams.globalAxis(:,Z),  'numberCheck', 1, 'returnType', true), ...
      HT_Face_Find(lMod_OvfSW_faces,            'orientation', -lHiveParams.globalAxis(:,Z), 'numberCheck', 1, 'returnType', true), ...
      struct(),...
      setfield(lOptions, 'strictMatch', false));
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);
  endfor


  % ========================== Radiation inside hive body ========================
  % At this stage, all faces around the hive body
  if lRadiation.enable_hiveBody
    lMod_HiveBody_intFaces = cell(6,1);

    for i=1:4
      lMod_HiveBody_intFaces{i} = HT_Face_Union(...
                                      sprintf('mSW1_%s', FaceDIRstrL{i}), ...
                                      {lMod_SideWall1_facesCell{i}(FaceZM), lMod_SideWall2_facesCell{i}(FaceZM)}, ...
                                      'coplanar', true, 'allowIntersect', false, 'mesh', true);
    endfor

    lMod_HiveBody_intFaces{5} = lMod_SideWallBot_faces{FaceZP};
    lMod_HiveBody_intFaces{6} = lMod_Overframe_faces{FaceZM};

    lMod_Tmp = HT_Model_RadiationInBox(...
                  "mHB_rad", ...
                  struct( 'face', { lMod_HiveBody_intFaces }, ...
                          'emissivity', [ HT_Material_GetEmissivity(lHiveMat.sidewalls_int, HT_WAVELENGTH_IR_MED) * ones(5,1); ...
                                          HT_Material_GetEmissivity(lHiveMat.overframe_int, HT_WAVELENGTH_IR_MED)]), ...
                  lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);
  endif

  % ============================= Underroof layer =================================
  [lMod_Underroof, lMod_Underroof_faces] = HT_Model_Conduction3D(...
    'mUnderroof', ...
    struct( 'base',     lMod_Overframe_faces{FaceZP},     ...
            'length',   lHiveParams.underroofThickness,   ...
            'material', lHiveMat.underroof,               ... % Retrieve the material specified in the header (generally air)
            'n',        [NA NA lMesh.overframe(Z)],       ...
            'zGridType', 'ff'...
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Underroof, lOptions);

  if lPlotOptions.plotGeometry
    HT_Plot_Face(lMod_Underroof_faces, l3DPlotOptions);
  endif

  % Reorder the face to match global axis
  % lMod_Underroof_faces is now a cell array instead of a struct array
  lMod_Underroof_faces = HT_Face_Find(lMod_Underroof_faces, 'orientation', FaceDIRmat, 'numberCheck', 1, 'returnType', true);

  % Connect underroof to the overframe
  lMod_Tmp = HT_Model_ConnectFaces(...
    'mcOverframeUnderroof', ...
    lMod_Overframe_faces{FaceZP}, ...
    lMod_Underroof_faces{FaceZM},...
    struct(),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);


  lMod_UnderSideWall_facesCell = cell(4,1);    % Backup faces for later use
  % Add surrounding wood sides
  for i=1:4
    % Retrieve the base face of overframe layer

    [lMod_UnderSW, lMod_UnderSW_faces] = HT_Model_Conduction3D(...
      sprintf('mUnderSW_%s', FaceDIRstrL{i}),         ... % Model name
      struct( 'base',     lMod_Underroof_faces{i},    ... % Base face over which the 1D model is built
              'length',   lHiveParams.wallThickness,  ...
              'material', lHiveMat.sidewalls,         ...
              'n',        [NA NA lMesh.sidewalls],    ... % NA is used since mesh is extracted from face <base>
              'zGridType', 'hh',                      ... % Half nodes are used on boundaries
              'xGridType', 'ff',...
              'yGridType', 'ff'...
              ),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_UnderSW, lOptions);
    lMod_UnderSideWall_facesCell{i} = lMod_UnderSW_faces;

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_UnderSW_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Connect underroof sidewalls to the underroof
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcUnderSW_%s', FaceDIRstrL{i}),  ...
      lMod_Underroof_faces{i},                  ...
      lMod_UnderSW_faces(FaceZM),               ... % For 1D conduction problem, FaceXM is always the base face
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect underroof sidewalls to overframe sidewalls
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcUnderSW_OvfSW_%s', FaceDIRstrL{i}),  ...
      HT_Face_Find(lMod_OvfSideWall_facesCell{i}, 'orientation', lHiveParams.globalAxis(:,Z),  'numberCheck', 1, 'returnType', true), ...
      HT_Face_Find(lMod_UnderSW_faces,            'orientation', -lHiveParams.globalAxis(:,Z), 'numberCheck', 1, 'returnType', true), ...
      struct(),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);
  endfor


  % =================================== Roof =====================================
  lRoofSize = lHiveParams.dimension + 2*lHiveParams.wallThickness + 2*lHiveParams.roofIntExtra;
  lRoofPosition = lHiveParams.globalPosition + ...
                  [ -lHiveParams.wallThickness-lHiveParams.roofIntExtra(X);...
                    -lHiveParams.wallThickness-lHiveParams.roofIntExtra(Y);...
                    lMod_Underroof_faces{FaceZP}.globalPosition(Z)];


  [lMod_Roof, lMod_Roof_faces] = HT_Model_Conduction3D(...
    'mRoof', ...
    struct( 'globalPosition', lRoofPosition,                    ...
            'dim',            lRoofSize,                        ...
            'axis',           lHiveParams.globalAxis,                      ...
            'length',         lHiveParams.roofThickness,        ...
            'material',       lHiveMat.roof_ext,                ... % Retrieve the material specified in the header
            'n',              [NA NA lMesh.roof],                        ... % Only specify the number of nodes for in-depth direction
            'zGridType',      'ff',                             ... % Specify full node size on boundaries (instead of half). This mode must be selected to allow 1 node only in the z-direction
            'xGrid',          { lMod_Underroof_faces{FaceZP} }, ... % Underroof faces are provided as "xgrid" so that the xdirection mesh will be adpated
            'yGrid',          { lMod_Underroof_faces{FaceZP} }, ... % Underroof faces are provided as "xgrid" so that the xdirection mesh will be adpated
            'xGridType',      'ff',                             ...
            'yGridType',      'ff'
            ),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Roof, lOptions);

  lMod_Tmp = HT_Model_Connect('mRoof_ExtAir', lAirExtNode,                  ... % Node of outside air
                                    lMod_Roof_faces(FaceZP),                ... % Top face of the roof
                                    'g', lConvection.roof_ext{FaceZP},      ... % Convection coefficient specified by user
                                    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Connect underroof to roof
  lMod_Tmp = HT_Model_ConnectFaces(...
    'mcUnderRoof_Roof',  ...
    lMod_Underroof_faces{FaceZP}, ...
    lMod_Roof_faces(FaceZM), ...
    struct('intersect', true),...
    lOptions);
  lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  % Plot geometry
  if lPlotOptions.plotGeometry
    HT_Plot_Face(lMod_Roof_faces, setfield(l3DPlotOptions, 'explodeOffset', [0 0 0.1])); %, l3DPlotOptions); % Display the geometry
  endif


  % Reorder the face to match global axis
  % lMod_Underroof_faces is now a cell array instead of a struct array
  lMod_Roof_faces = HT_Face_Find(lMod_Roof_faces, 'orientation', FaceDIRmat, 'numberCheck', 1, 'returnType', true);

  % =============================== Roof sides ===================================
  lRoofSizeAxis = { FaceDIRmat(:, [FaceYP,FaceZM,FaceXM]), ...
                    FaceDIRmat(:, [FaceYM,FaceZM,FaceXP]), ...
                    FaceDIRmat(:, [FaceXM,FaceZM,FaceYM]), ...
                    FaceDIRmat(:, [FaceXP,FaceZM,FaceYP]) };
  lRoofSide_faces = cell(4,1);    % Store roof side faces for future use
  lAirSide_faces = cell(4,1);

  for i=1:4
    [lMod_RoofSide, lMod_RoofSide_faces] = HT_Model_Conduction3D(...
      sprintf('mRoofSide_%s', FaceDIRstrL{i}),  ...
      struct( 'base',           lMod_Roof_faces{i},               ... % Base face
              'globalPosition', [NA; NA; lRoofPosition(Z)],       ...
              'dim',            [NA; lHiveParams.roofHeight],     ... %
              'axis',           lRoofSizeAxis{i},                 ...
              'length',         lHiveParams.roofThickness,        ...
              'material',       lHiveMat.roofsides_ext,           ... % Retrieve the material specified in the header
              'n',              [NA NA lMesh.roof],               ...
              'zGridType',      'ff',                             ... % Specify full node size on boundaries (instead of half). This mode must be selected to allow 1 node only in the z-direction
              'xGrid',          [],                               ... % Load from <.base> face
              'xGridType',      'ff',                             ... % No half node on boundaries
              'yGrid',          { HT_Face_RemoveMesh({ lMod_Underroof_faces{i},     ... % Underroof, overframe, and hive body
                                                       lMod_Overframe_faces{i},     ... % are used to split the roof side along
                                                       lMod_HiveBody2_faces(i)}) }, ... % the z-direction
              'yGridType',      'ff' ...
              ),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_RoofSide, lOptions);
    lRoofSide_faces{i} = lMod_RoofSide_faces;

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_RoofSide_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Connect the roof side to the roof
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcRoofSide_Roof_%s', FaceDIRstrL{i}),      ...
      lMod_RoofSide_faces(FaceYM),                        ...
      lMod_Roof_faces{i},                                 ...
      struct( 'uv1axis', lRoofSizeAxis{i}(:,[X Z]),       ...
              'uv2axis', lRoofSizeAxis{i}(:,[X Y])),      ...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Build the air layer between the roof side and the hive
    % The union of faces (underroof, overframe, hive body 2, hive body 1)
    % is done.
    lSideFace = HT_Face_Union(sprintf('sideUnion_%s', FaceDIRstrL{i}), ...
                                   {  lMod_SideWall1_facesCell{i}(FaceZP),...
                                      lMod_SideWall2_facesCell{i}(FaceZP),...
                                      lMod_OvfSideWall_facesCell{i}(FaceZP), ...
                                      lMod_UnderSideWall_facesCell{i}(FaceZP) }, ...
                                      'mesh', false,        ... % Setup the maximum restrictive conditions
                                      'allowIntersect', false,  ...
                                      'coplanar', true);

    % Compute the intersection of the previous union and the internal face of the roof side
    % Even if they are not coplanar, the function preserves the position of the first
    % face in the list.
    lSideFace = HT_Face_Intersect(sprintf('sideIntersect_%s', FaceDIRstrL{i}), ...
                                   { lSideFace, lRoofSide_faces{i}(FaceZM) }, ...
                                   'mesh', false, ...
                                   'coplanar', false); ... % Roof side face and hive faces are not coplanar

    [lMod_AirSide, lMod_AirSide_faces] = HT_Model_Conduction3D(...
      sprintf('mAirSide_%s', FaceDIRstrL{i}),  ...
      struct( 'base',           lSideFace,                                ... % Internal face
              'axis',           lRoofSizeAxis{i},                         ...
              'length',         lHiveParams.roofIntExtra([1 2 1 2])(i),   ...
              'material',       lHiveMat.roofsideAir,                     ... % Retrieve the material specified in the header
              'n',              [NA NA lMesh.airside],                    ...
              'zGridType',      'ff',                                     ... % Specify full node size on boundaries (instead of half). This mode must be selected to allow 1 node only in the z-direction
              'xGrid',          { lRoofSide_faces{i}(FaceZM) },           ... % Load from <.base> face
              'xGridType',      'ff',                                     ... % No half node on boundaries
              'yGrid',          { lRoofSide_faces{i}(FaceZM) },           ... % the z-direction
              'yGridType',      'ff' ...
              ),...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_AirSide, lOptions);
    lAirSide_faces{i} = lMod_AirSide_faces;

    % Plot geometry
    if lPlotOptions.plotGeometry
      HT_Plot_Face(lMod_AirSide_faces, l3DPlotOptions); %, l3DPlotOptions); % Display the geometry
    endif

    % Connect the side air layer to the roof side
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcRoofSide_AirSide_%s', FaceDIRstrL{i}),   ...
      lMod_AirSide_faces(FaceZP),                         ...
      lRoofSide_faces{i}(FaceZM),                         ...
      struct('intersect', true),                          ...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect the side air layer to the roof
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcRoof_AirSide_%s', FaceDIRstrL{i}),               ...
      lMod_AirSide_faces(FaceYM),                                 ...
      lMod_Roof_faces{FaceZM},                                    ...
      struct('intersect', true),                                  ...
      lOptions);
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect the side air layer to the underroof side walls
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcAirSide_UnderSW_%s', FaceDIRstrL{i}),              ...
      lMod_AirSide_faces(FaceZM),                                   ...
      lMod_UnderSideWall_facesCell{i}(FaceZP),                      ...
      struct('intersect', true),              ...
      setfield(lOptions, 'strictMatch', false));
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect the side air layer to the overframe side walls
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcAirSide_OvfSW_%s', FaceDIRstrL{i}),                ...
      lMod_AirSide_faces(FaceZM),                                   ...
      lMod_OvfSideWall_facesCell{i}(FaceZP),                        ...
      struct('intersect', true),              ...
      setfield(lOptions, 'strictMatch', false));
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect the side air layer to the hive body 2
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcAirSide_SW2_%s', FaceDIRstrL{i}),                  ...
      lMod_AirSide_faces(FaceZM),                                   ...
      lMod_SideWall2_facesCell{i}(FaceZP),                          ...
      struct('intersect', true),              ...
      setfield(lOptions, 'strictMatch', false));
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

    % Connect the side air layer to the hive body 1
    lMod_Tmp = HT_Model_ConnectFaces(...
      sprintf('mcAirSide_SW1_%s', FaceDIRstrL{i}),                  ...
      lMod_AirSide_faces(FaceZM),                                   ...
      lMod_SideWall1_facesCell{i}(FaceZP),                          ...
      struct('intersect', true),              ...
      setfield(lOptions, 'strictMatch', false));
    lModel = HT_Model_Merge(lModel, lMod_Tmp, lOptions);

  endfor

  clear lRoofSize lRoofPosition;

  disp('==============================================');
  disp(sprintf('Model built in %.2f s', lModel.timer.buildingTime));

  M = lModel;
  faces = struct();
  faces.hiveBody1 = lMod_HiveBody1_faces;
  faces.hiveBody2 = lMod_HiveBody2_faces;
  faces.roofSides = lRoofSide_faces;
  faces.roof = lMod_Roof_faces;
  faces.sideWall1 = lMod_SideWall1_facesCell;
  faces.sideWall2 = lMod_SideWall2_facesCell;

  volumes = struct();
  volumes.hiveBody = lMod_HiveBodyVolumes;

endfunction
