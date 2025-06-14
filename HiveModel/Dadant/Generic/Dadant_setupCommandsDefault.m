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

function [lCmd lCmdInfos]= Dadant_setupCommandsDefault(lParams, lFaces, lOptions)
  HT_ImportConstants();

  % ========================================================================
  %                          PRELIMINAR COMPUTATIONS
  % ========================================================================
  % Retrieve the sun position at current date time
  [zenithAngle azimuthAngle] = HT_GetSolarAngle(...
                                  lParams.hiveParams.date(1),                    ... Year
                                  lParams.hiveParams.date(2:3),                  ... [Month Day]
                                  sum(lParams.hiveParams.time .* [3600 60 1]) +  ...
                                    lParams.computation.startTime + (0:(lParams.computation.nt-1)) * lParams.computation.timeStep,          ... Day time
                                  lParams.hiveParams.longitude,                  ... Longitude
                                  lParams.hiveParams.latitude);                  %  Latitude

  % Retrieve the cosine of angle between each hive faces and the sun
  % Each row of < cosBeta (dim Nx6) > refers to a different time step
  % XM, XP, YM, YP, ZM, ZP (top)
  [cosBeta] = HT_GetHiveSolarAngles(...
                            zenithAngle,                                  ... % Sun zenith angle (90°=sun at horizon)
                            azimuthAngle,                                 ... % Sun azimuth angle (radian)
                            [lParams.hiveParams.azimuth, lParams.hiveParams.tilt],      ...
                            struct('angle', true));                       % angle=true means previous
                                                                          %  argument is expressed as [azimuth tilt]

  % Roof and roof sides: Retrieve the exposition coefficient (negative values are clamped to 0)
  lRoofExpositionMat = HT_GetSolarFluxFactor(   ...
                            cosBeta,            ... % Cosine between face and sun for each faces
                            [HT_Material_GetEmissivity(lParams.hiveMat.roofsides_ext, HT_WAVELENGTH_VISIBLE)*ones(1,4), ... XM, XP, YM, YP
                              0, HT_Material_GetEmissivity(lParams.hiveMat.roof_ext, HT_WAVELENGTH_VISIBLE)]);

  % Side walls: Exposition coefficient
  lHiveBodySWExpositionMat = HT_GetSolarFluxFactor(   ...
                            cosBeta,                  ... % Cosine between face and sun for each faces
                            [HT_Material_GetEmissivity(lParams.hiveMat.sidewalls, HT_WAVELENGTH_VISIBLE)*ones(1,4), ... XM, XP, YM, YP
                              0, 0]);

  % Nullify exposition coefficients when the sun is below the horizon
  % i.e. when zenith angle is higher than 90°=pi/2
  ##lExpositionMat(zenithAngle >= pi/2, :) = 0;
  lZenithAttenuation = 0.8; % Voir Solar Irradiance on Wikipedia (radiation at ground level over/extra atmosphere) with at zenith
  lSolarIrradiance = 1361 * exp(log(lZenithAttenuation) ./ max(0.0, cos(zenithAngle)));
  % Multiply each line by the solar irradiance

  lHiveBodySWExpositionMat = lSolarIrradiance .* lHiveBodySWExpositionMat;
  lRoofExpositionMat = lSolarIrradiance .* lRoofExpositionMat;

  lCmdInfos.roofIrradiance = lRoofExpositionMat;
  lCmdInfos.hiveBodyIrradiance = lHiveBodySWExpositionMat;

  clear lAreaFaces lSolarIrradiance lZenithAttenuation cosBeta zenithAngle azimuthAngle; % lExpositionMat

  % ========================================================================
  %                              COMMANDS
  % ========================================================================
  lCmd = HT_Cmd_Init(
          'name', 'cmd_airext',           ...
            'type', 'temperature',        ...
            'nodes', 'nAirExt',           ...
            'data', 10,                   ...
            'enable', true);                ...

  % Add the roof sides (solar irradiance)
  for i=1:4
    lCmd = [lCmd; HT_Cmd_Init('name', sprintf('cmd_irradiance_roof_%s', FaceDIRstrL{i}), ...
                              'type', 'flux', ...
                              'nodes', { lFaces.roofSides{i}(FaceZP) } , ...
                              'data', lRoofExpositionMat(:,i)', ...
                              'enable', true)];
  endfor

  % Add the roof top (solar irradiance)
    lCmd = [lCmd; HT_Cmd_Init('name', sprintf('cmd_irradiance_roof_%s', FaceDIRstrL{FaceZP}), ...
                              'type', 'flux', ...
                              'nodes', { lFaces.roof{FaceZP} } , ...
                              'data', lRoofExpositionMat(:,FaceZP)', ...
                              'enable', true)];

  % Add the hive sides sides (solar irradiance)
  for i=1:4
    lCmd = [lCmd; HT_Cmd_Init('name', sprintf('cmd_irradiance_sw_%s', FaceDIRstrL{i}), ...
                              'type', 'flux', ...
                              'nodes', { lFaces.sideWall1{i}(FaceZP) } , ...
                              'data', lHiveBodySWExpositionMat(:,i)', ...
                              'enable', true)];
  endfor

endfunction
