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

% Retrieve the cosine of angle between each 8 hive faces and the sun
% Each row of < cosBeta (dim Nx8) > refers to a different time step
lNormalVectorList = cell2mat(cellfun(@(v) v.norm', lLateralOutsideFaces, 'UniformOutput', false));
lNormalVectorList = [lNormalVectorList; [0 0 1]]; % The vector towards Z+ for the roof

IRRADIANCE_COLUMN_ROOF = 9;

[cosBeta] = HT_GetHiveSolarAngles(...
                          zenithAngle,                                  ... % Sun zenith angle (90°=sun at horizon)
                          azimuthAngle,                                 ... % Sun azimuth angle (radian)
                          lNormalVectorList, ...
                          struct('angle', false));                       % angle=true means previous
                                                                        %  argument is expressed as [azimuth tilt]

lZenithAttenuation = 0.8; % Voir Solar Irradiance on Wikipedia (radiation at ground level over/extra atmosphere) with at zenith
lSolarIrradiance = 1361 * exp(log(lZenithAttenuation) ./ max(0.0, cos(zenithAngle'))); % Row vector

lSunAbsorptivityVec = [ repmat(HT_Material_GetEmissivity(lHiveMat.extWalls, HT_WAVELENGTH_VISIBLE), 1, 8) ...
                        HT_Material_GetEmissivity(lHiveMat.roof, HT_WAVELENGTH_VISIBLE) ];
lSunIrradianceMatrix = (max(cosBeta,0) .* lSunAbsorptivityVec) .* lSolarIrradiance';

% This matrix will contain the sun direction (in column) at each time step
##lSunDirectionMatrix = [cos(azimuthAngle'+pi); sin(azimuthAngle'+pi); cos(zenithAngle')];
##lSunDirectionMatrix(1:2,:) = lSunDirectionMatrix(1:2,:) .* sin(zenithAngle');

% Split the matrix J [nodes x angles] into parts to match each sides of the cylinder (lMod_CylinderFaces.outside)
% Each cell will contain <lMesh.bodyHeight> nodes.
##Jcell = mat2cell(J.*lSolarIrradiance, repmat(lMesh.bodyHeight, 1, lMesh.circular));


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
          'nodes',  [lRoofOutsideFaces; lMod_AirCylinderFaces(FaceXP)], ...
          'data', lSunIrradianceMatrix(:,IRRADIANCE_COLUMN_ROOF)');

for i=1:8
  lCmd = [lCmd; HT_Cmd_Init(...
                  'name', sprintf('cmd_side_%d', i), ...
                  'type', 'flux', ...
                  'nodes', lLateralOutsideFaces(i), ...
                  'data', lSunIrradianceMatrix(:,i)')];
endfor
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


