
%*****************************************************************
% Title: sun_vector
% Author: Lisa Reeh
% Date written: 14 Feb. 2002
% Date last modified: 27 May 2003
% Copyright 2002 Univesity of Colorado, Boulder
%
% Purpose: computes the vector from the center of the earth
% to the center of the sun
%
% Input:
%       T: number of elapsed Julian centuries (J2000 epoch)
%
% Output:
%       r_sun_vector: ECI sun vector (meters)
%
% Reference: Vallado, p 183 (1st ed.), OR p.263 (2nd ed.)
%*****************************************************************

function r_sun_vector = sun_vector(T)

d2r = pi/180;
AU = 149597870700;  % Astronomical Unit [in meters]

% determine the mean longitude of the sun in degrees
lambda = mod( (280.4606184 + 36000.77005361 * T), 360 );

% determine the mean anomaly of the sun in degrees
M = mod( (357.5277233 + 35999.05034 * T), 360 );

% Calculate the ecliptic longitude in degrees
lambda_ecl = mod(lambda + 1.914666471 * sin(M*d2r) + 0.019994643 * sin(2*M*d2r), 360);

% compute the magnitude of the sun vector
r_sun_mag = 1.000140612 - 0.016708617*cos(M*d2r) - 0.000139589*cos(2*M*d2r);

% approximate the obliquity of the ecliptic
ob_ecl = mod(23.439291 - 0.0130042 * T, 360);

% Compute sun vector in ECI frame
rm = repmat(r_sun_mag,1,3);  %use this to remove the requirement of a loop

% Lots of transposes to keep matrix dimensions consistant
r_sun_vector = AU*rm.*([cos(lambda_ecl.*d2r)'; cos(ob_ecl.*d2r)'.*sin(lambda_ecl.*d2r)';...
    sin(ob_ecl.*d2r)'.*sin(lambda_ecl.*d2r)'])';

end