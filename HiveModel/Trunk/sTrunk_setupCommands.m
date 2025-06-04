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
%                          PRELIMINAR COMPUTATIONS
% ========================================================================
timeVector = sum(lHiveParams.time .* [3600 60 1]) + (0:(lComputation.nt-1))' * lComputation.timeStep;

% ==============================================================================
% Retrieve the sun position at current date time
% Return angles are in Radian unless option 'unit' is specified.
% azimuthAngle is relative to north and ccw (-90°=east)
% zenithAngle is relative to vertical
[zenithAngle azimuthAngle] = HT_GetSolarAngle(...
                                lHiveParams.date(1),                    ... Year
                                lHiveParams.date(2:3),                  ... [Month Day]
                                timeVector,                             ... Day time
                                lHiveParams.longitude,                  ... Longitude
                                lHiveParams.latitude);                  %  Latitude

lZenithAttenuation = 0.8; % Voir Solar Irradiance on Wikipedia (radiation at ground level over/extra atmosphere) with at zenith
lSolarIrradiance = 1361 * exp(log(lZenithAttenuation) ./ max(0.0, cos(zenithAngle'))); % Row vector

% This matrix will contain the sun direction (in column) at each time step
lSunDirectionMatrix = [cos(azimuthAngle'+pi); sin(azimuthAngle'+pi); cos(zenithAngle')];
lSunDirectionMatrix(1:2,:) = lSunDirectionMatrix(1:2,:) .* sin(zenithAngle');

% ==============================================================================
% Compute the projected shadow of the roof on the hive main part (cylinder)

% Retrieve the node positions along the vertical axis of the hive
lWallFace = HT_Face_AdjustAxis(lMod_CylinderFaces.outside(1), lHiveParams.globalAxis(:,Z));
lWallVertices = HT_Face_GetAbsoluteData(lWallFace, "xset");

% Azimuth angle passed to this function is relative to the first basis vector
% Passing azimuth=0 means the sun is oriented to south.
% Passing azimuth=pi/2 means the sun is oriented to east
[~, J] = HT_Shadow_GetDiskOverCylinder(...
                    lHiveParams.stoneRadius,                                ... % Radius of the roof
                    lHiveParams.bodyRadius + lHiveParams.wallThickness,     ...
                    lWallVertices,                                          ...
                    (2*pi/lMesh.circular)*((0:lMesh.circular)' - 0.5),      ... % Copy the way the disk is subdivided by Model_ConductionCylinder3D
                    [azimuthAngle'+pi ; pi/2-zenithAngle'],                   ... % Azimuth and elevation angle. Pi shift the reference of the azimuth, and pi/2 converts zenith angle to elevation angle
                      'deltaDisk', lHiveParams.roofShift, ...
                      'epsType', lHiveMat.walls.epsModel, ...
                      'power', lHiveMat.walls.epsParameters, ...
                      'sort', 'height');                                    ... % We ask the function to put values in matrix J based on height first (then theta)
                                                                                % since outside faces are naturally ordered by height

clear lWallFace;

% Split the matrix J [nodes x angles] into parts to match each sides of the cylinder (lMod_CylinderFaces.outside)
% Each cell will contain <lMesh.bodyHeight> nodes.
Jcell = mat2cell(J.*lSolarIrradiance, repmat(lMesh.bodyHeight, 1, lMesh.circular));


##% Roof and roof sides: Retrieve the exposition coefficient (negative values are clamped to 0)
##lRoofExpositionMat = HT_GetSolarFluxFactor(   ...
##                          cosBeta,            ... % Cosine between face and sun for each faces
##                          [HT_Material_GetEmissivity(lHiveMat.roofsides_ext, HT_WAVELENGTH_VISIBLE)*ones(1,4), ... XM, XP, YM, YP
##                            0, HT_Material_GetEmissivity(lHiveMat.roof_ext, HT_WAVELENGTH_VISIBLE)]);
##
##% Side walls: Exposition coefficient
##lHiveBodySWExpositionMat = HT_GetSolarFluxFactor(   ...
##                          cosBeta,                  ... % Cosine between face and sun for each faces
##                          [HT_Material_GetEmissivity(lHiveMat.sidewalls, HT_WAVELENGTH_VISIBLE)*ones(1,4), ... XM, XP, YM, YP
##                            0, 0]);
##
##% Nullify exposition coefficients when the sun is below the horizon
##% i.e. when zenith angle is higher than 90°=pi/2
##lExpositionMat(zenithAngle >= pi/2, :) = 0;
##% Multiply each line by the solar irradiance
##lExpositionMat = lSolarIrradiance .* lExpositionMat;
##% Retrieve the area of each external faces
####lAreaFaces = HT_Face_GetAbsoluteData(cellfun(@(v) v(FaceXP), lWallFaces, 'UniformOutput', false), 'area');
##% Multiply each columns by the face area
####lSolarRadiationMat = lExpositionMat .* cell2mat(lAreaFaces)';
##
##clear lAreaFaces lExpositionMat lSolarIrradiance lZenithAttenuation zenithAngle azimuthAngle;

% ========================================================================
%                              COMMANDS
% ========================================================================
lCmd = HT_Cmd_Init(
        'name', 'cmd_airext', ...
          'type', 'temperature', ...
          'nodes', 'nAirExt', ...
          'data', 10, ...
##        'name', 'cmd_sun_wall', ...
##          'type', 'flux', ...
##          'nodes', lMod_CylinderFaces.outside, ...
##          'data', Jcell, ...
##          'area', 'none',...
        'name', 'cmd_sun_roof', ...
          'type', 'flux',...
          'nodes',  lMod_RoofStoneFaces.top, ...
          'data', 1000); %repmat(1000, numel(lMod_RoofStoneFaces.top), 1)); %lSolarIrradiance); % .* max(0, lMod_RoofStoneFaces.top.norm' * lSunDirectionMatrix)); %, ... % dot product between normal and sun direction
##          'data', repmat(1000, 1, numel(timeVector))); %lSolarIrradiance); % .* max(0, lMod_RoofStoneFaces.top.norm' * lSunDirectionMatrix)); %, ... % dot product between normal and sun direction
##        'name', 'cmd_sun_roofsides', ...
##          'type', 'flux', ...
##          'nodes', lMod_RoofStoneFaces.outside, ...
##          'data', arrayfun(@(f) lSolarIrradiance .* max(0, f.norm' * lSunDirectionMatrix), lMod_RoofStoneFaces.outside, 'uniformoutput', false)); % lSolarRadiationMat(:,FaceXM), ...
##        'name', 'cmd_sun_xp', ...
##          'type', 'flux', ...
##          'nodes', lRoofSide_faces{FaceXP}(FaceZP).nodes, ...
##          'data', lSolarRadiationMat(:,FaceXP), ...
##        'name', 'cmd_sun_ym', ...
##          'type', 'flux', ...
##          'nodes', lRoofSide_faces{FaceYM}(FaceZP).nodes, ...
##          'data', lSolarRadiationMat(:,FaceYM), ...
##        'name', 'cmd_sun_yp', ...
##          'type', 'flux', ...
##          'nodes', lRoofSide_faces{FaceYP}(FaceZP).nodes, ...
##          'data', lSolarRadiationMat(:,FaceYP), ...
##        'name', 'cmd_sun_zp',     ...
##          'type', 'flux',         ...
##          'nodes', lMod_Roof_faces(FaceZP).nodes, ... Roof face
##          'data', lSolarRadiationMat(:,FaceZP));


