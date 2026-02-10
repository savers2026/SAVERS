%****************************************************************************
% Title: compute_shadow
% Author: Lisa Reeh
% Date written: Jan. 2003
% Date last modified: 27 May 2003
% Copyright 2002 University of Colorado, Boulder
%
% Purpose: Determines if the current GPS satellite is in the earth's
%       penumbra region.
%
% Input: 
%       r_sun_vector: vector from earth to sun (m, ECI)
%       sat_pos_vect: position vector of the satellite (m, ECI)
%
% Output:
%       shadow: (1) if satellite is in earth penumbra
%               (0) if satellite is not in earth penumbra
%
% Reference: Vallado, 2nd ed., Algorithm 34, with slight modifications
%****************************************************************************

function shadow = compute_shadow(r_sun_vector, sat_pos_vect)

R_sun = 695990000;      % sun's radius, m (Conway & Prussing)    
R_earth = 6378136.3;    % earth's radius, m (Vallado, 2nd ed., JGM-2)
r_sun_mag = norm(r_sun_vector);
sat_pos_mag = norm(sat_pos_vect);
angle_penumbra = atan2((R_sun + R_earth),r_sun_mag);

if dot(r_sun_vector, sat_pos_vect) < 0
    % determine the angle between r_sun_vector and sat_pos_vect
    angle = acos(dot(r_sun_vector, sat_pos_vect)/(r_sun_mag * sat_pos_mag));
    
    sat_horiz = sat_pos_mag * cos(angle);
    sat_vert = sat_pos_mag * sin(angle);
    pen_vert = R_earth + tan(angle_penumbra)*sat_horiz;
    
    if sat_vert <= pen_vert
        shadow = 1;
    else
        shadow = 0;
    end
else
    shadow = 0;
end
end