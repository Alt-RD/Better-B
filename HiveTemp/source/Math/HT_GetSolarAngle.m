%  This file is part of project HiveTemp.
%
%  Copyright (c) 2022 Montpellier-University, AltRD-Emmanuel Ruffio
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
% This function returns the position of the sun in the sky seen by an
% observer at ground level.
%
% Input arguments:
% year : the current year (ex: 2019)
% day : the day of the year [1; 365 (or 366 for leap years)] or dim 2=[month, day]  (month and day minimum values is 1)
% dayTime [dim Nx1]: utc time in second (can be negative)
% latitude: latitude of current position (real value [in deg] or dim 3 vector=[deg minute second])
% longitude: longitude of current position (real value [in deg] or dim 3 vector=[deg minute second])
%   Ex: Montpellier: Lat/Long: 43°37'00"N / 3°54'00"E   ou  (43.617° / 3.9°)
% options: struct('paramname', paramvalue,...) containing options
%   - 'unit' -> = none/'degres' specify the unit of output angles. If not specified, returned values are expressed in radian.
%
% Output arguments:
% -> phiDeg: [column vector] phiDeg is the zenith angle (0 means the solar is at zenith, and 90° means it is at the horizon level
%            [in radian, unless options.unit='degres']
% -> thetaDeg: [column vector] solar azimuth (relative to geo. North and counter-clockwise) [in radian, unless options.unit='degres']
% -> declDeg: [column vector] the declination of the sun. The solar declination varies
%             from -23.44° at the (northern hemisphere) winter solstice, through 0° at the
%             vernal equinox, to +23.44° at the summer solstice/
%             Angle between sun and equator.
function [phi theta decl] = HT_GetSolarAngle(year, day, dayTime, latitude, longitude, options)

if nargin < 6 || isempty(options), options = struct(); endif;

if ~isfield(options, "unit")
  options.unit = 0;
elseif strcmpi(options.unit, "degres")
  options.unit = 1;
else
  options.unit = 0;
endif

assert(any(numel(latitude) == [3 1]), 'Latitude must be a vector of dim 3 or 1');
assert(any(numel(longitude) == [3 1]), 'Longitude must be a vector of dim 3 or 1');

dayTime = dayTime(:); % Make sure the vector is a column vector

yearInt = int32(year);
lIsDiv4 = ((yearInt/4)*4 == yearInt);
lIsDiv100 = ((yearInt/100)*100 == yearInt);
lIsDiv400 = ((yearInt/400)*400 == yearInt);
lIsLeapYear = (lIsDiv4 && !lIsDiv100) || lIsDiv400;

% Convert [month day] to day index [0; (day of the year)
if columns(day) > 1
  assert(all(all(day > 0)), 'Months and days indices minimum is 1');
  dayMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
  if lIsLeapYear
    dayMonth(2) += 1; % Leap year
  endif

  dayVec = zeros(rows(day), 1);
  for i=1:numel(dayVec)
    dayVec(i) = sum(dayMonth(1:(day(i,1)-1))) + day(i, 2);
  endfor
  day = dayVec;
  clear dayVec, dayMonth;
else
  assert(all(day > 0), 'Day minimum value is 1');
endif

if numel(longitude) != 1 % Conversion des degrés°minutes en degré
  longitude = longitude(1) + longitude(2)/60 + longitude(3)/3600;
endif

if numel(latitude) != 1 % Conversion des dégrés°minutes en degré
  latitude = latitude(1) + latitude(2)/60 + latitude(3)/3600;
endif

% Convert latitude and longitude to radian
latRad = latitude/360*2*pi;
longRad = longitude/360*2*pi;

% Voir solareqns.PDF, General Position Position Calculations
% NOAA Global Monitoring Division

% First, the fractional year (?) is calculated, in radians
% It is obviously UTC hour since declinaison does not depends on our position.
if lIsLeapYear
	lDayCount = 366;
else
	lDayCount = 365;
endif

gamma = 2*pi/lDayCount * (day-1 + (dayTime/3600 - 12)/24);

% From ?, we can estimate the equation of time (in minutes)
eqtime = 229.18 * ( 0.000075 ...
                    +0.001868 * cos(gamma) ...
                    -0.032077 * sin(gamma) ...
                    -0.014615 * cos(2*gamma) ...
                    -0.040849 * sin(2*gamma) );

% the solar declination angle (in radians).
decl = 0.006918 -0.399912 * cos(gamma) ...
                +0.070257 * sin(gamma) ...
                -0.006758 * cos(2*gamma) ...
                +0.000907 * sin(2*gamma) ...
                -0.002697 * cos(3*gamma) ...
                +0.00148 * sin(3*gamma);
% https://www.esrl.noaa.gov/gmd/grad/solcalc/glossary.html
% solar declination - the declination of the sun. The solar declination varies
% from -23.44° at the (northern hemisphere) winter solstice, through 0° at the
% vernal equinox, to +23.44° at the summer solstice. The variation in solar
% declination is the astronomical description of the sun going south (in the
% northern hemisphere) for the winter. Click on Solar Declination Graph to see
% how the solar declination varies over the year. See Solar Paths Figure to see
% the seasonal solar paths projected on the celestial sphere. For a ground-based
% view of the seasonal solar paths for different latitudes, see: 0° (the Equator),
% 23°N (the Tropic of Cancer), 40°N (Boulder, CO), 71°N (the Arctic Circle),
% and 90° (the North Pole).

% Next, the true solar time is calculated in the following two equations.
% First the time offset is found, in minutes
time_offset = eqtime + 12/180*60*longitude;
% where eqtime is in minutes,
%       longitude is in degrees (positive to the east of the Prime Meridian),
%       timezone is in hours from UTC (U.S. Mountain Standard Time = –7 hours).

% and then the true solar time, in minutes.
% tst = hr*60 + mn + sc/60 + time_offset; % where hr is the hour (0 - 23), mn is the minute (0 - 59), sc is the second (0 - 59).
tst = dayTime/60 + time_offset;

% The solar hour angle, in radians, is:
ha = mod((tst/4 - 180) / 360 * 2*pi + pi, 2*pi) - pi; %radian

% The solar zenith angle (phi) can then be found from the hour angle (ha),
% latitude (lat) and solar declination (decl) using the following equation:
cosPhi = sin(latRad)*sin(decl) + cos(latRad)*cos(decl).*cos(ha);
phi = acos(cosPhi);

%[tst ha phiDeg]

% And the solar azimuth (?, degrees clockwise from north) is found from:
%cos180theta = -(sin(latitude) * cos(phi) - sin(decl)) ./ (cos(latitude)*sin(phi));

theta = zeros(size(gamma));
lPos = ha > 0;
lNeg = ~lPos;
theta(lPos) = mod(acos((sin(latRad)*cosPhi(lPos)-sin(decl(lPos)))./(cos(latRad)*sin(phi(lPos))))+pi, 2*pi);
theta(lNeg) = mod(3*pi-acos(((sin(latRad)*cosPhi(lNeg))-sin(decl(lNeg)))./(cos(latRad)*sin(phi(lNeg)))), 2*pi);
##theta(lPos) = acos((sin(latRad)*cosPhi(lPos)-sin(decl(lPos)))./(cos(latRad)*sin(phi(lPos))))+pi;
##theta(lNeg) = 3*pi-acos(((sin(latRad)*cosPhi(lNeg))-sin(decl(lNeg)))./(cos(latRad)*sin(phi(lNeg))));
theta = -theta;

if (options.unit == 1)
  decl = decl/(2*pi)*360;
  phi = phi/(2*pi)*360;
  theta = theta/(2*pi)*360;
endif

%if PhiDeg >0
%  theta = mod(acos((sin(latRad)*cosPhi-sin(decl))/(cos(latRad)*SIN(phi)))+pi;2*pi);
%else
%  theta = mod(3*pi-acos(((sin(latRad)*cos(decl))-sin(decl))/(cos(latRad)*SIN(phi))), 2*pi);
%endif
