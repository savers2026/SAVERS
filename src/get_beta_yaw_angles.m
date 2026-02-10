% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The Intellectual Property of the contents of this directory, code, data
% and files belong to SSTL and are @copyright SSTL 2018 unless otherwise
% specifically indicated within the code or file.
%
% The contents are being transferred to ESA under the agreement within
% Contract 4000122597//17/NL/CRS/hh - Design and development of a tool for 
% antenna carrier characterisation from flying space GNSS receiver
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Inner Functions: 
%       - sun_vector.m, to compute the Sun-pointing vector  
%
%  Inputs: 
%       - pos_vect: (3xn) position vector of the satellite (m, ECI)
%       - vel_vect: (3xn) velocity vector of the satellite (m/s, ECI)
%       - T: (nx1) number of elapsed Julian centuries (J2000 epoch)
%       - time: (nx1) measurement times (secs from epoch)
%       - tx_ECI_pos_unit: (3xn) unit vector of the satellite position (ECI)
%
%  Outputs: 
%       - yaw: the yaw angle of the satellite [deg] (see user manual for definition)
%       - alpha: Angle ALPHA [rad] (see user manual for definition)
%       - beta: % Angle BETA [deg] (see user manual for definition)
% 
%  References:
%       Original function information
%               Name: get_yaw_angle
%               Author: Lisa Reeh
%               Date written: 8 Feb. 2002
%               Date last modified: 27 May 2003
%               Copyright 2002 University of Colorado, Boulder
%       Yaw computation follows: 
%               Title: GPS Modeling and Analysis - Summary of Research: GPS Satellite Axial Ratio Predictions 
%               Author: Penina Axelrad and Lisa Reeh
%               Date: 23 Sept. 2002
%               Colorado Center for Astrodynamics Research - University of
%               Colorado, UCB 431
% ------------------------------------------------------------------------

function [yaw,alpha,beta] = get_beta_yaw_angles(pos_vect, vel_vect, T, time, tx_ECI_pos_unit)
yaw_rate = nan(length(time),1);
yaw = nan(length(time),1);
yaw_all = nan(length(time),1);

d2r = pi/180;                                               % convert from degrees to radians
r2d = 180/pi;                                               % convert from radians to degrees

%% BETA ANGLE
% Compute the Earth-Sun vector and their relative norm (in ECI frame)
r_sun_vector = sun_vector(T);                                   % Sun-Pointing vector    
r_sun_vector_norm = sqrt(sum(abs(r_sun_vector).^2,2));          % Magnitude of the Sun-Pointing vector (for each column)
r_sun_vector_norm = repmat(r_sun_vector_norm,1,3);
r_sun_unit = r_sun_vector./r_sun_vector_norm;                   % Sun-Pointing unit vector

% Compute the Orbit Norm vector (in ECI frame)
orb_normal = cross(pos_vect', vel_vect', 2);                    % Cross-Product: <GPS Position ,GPS Velocity> vectors (h vector)
orb_normal_norm = sqrt(sum(abs(orb_normal).^2,2));              % Magnitude of the Orbit Norm vecotr (for each column)
orb_normal_norm = repmat(orb_normal_norm,1,3);
orb_norm_unit = orb_normal./orb_normal_norm;                    % Orbit Norm unit vector

% Dot-Product <Orbit Norm unit, Sun-Pointing unit>
dot_prod = dot(orb_norm_unit,r_sun_unit,2);                     % Dot-Product: <Orbit Norm UNIT,Sun Pointing UNIT> vectors 

% Angle BETA [deg]
%beta = 90-rad2deg(acos(dot_prod));

% % Angle BETA [deg] NASA
orb_norm_unit_norm = vecnorm(orb_norm_unit,2,2);
r_sun_unit_norm = vecnorm(r_sun_unit,2,2);
beta = asin(dot_prod./(orb_norm_unit_norm.*r_sun_unit_norm))*r2d;

%% ALPHA ANGLE
% Vector "n", perpendicular to the projection of the sun vector onto the orbit plane and the satellite position vector
h_cross_rsun = cross(orb_normal,r_sun_vector,2);                % n-vector
h_cross_rsun_norm = sqrt(sum(abs(h_cross_rsun).^2,2));          % Magnitude of the n-vector (for each column)
h_cross_rsun_norm = repmat(h_cross_rsun_norm,1,3);              
h_cross_rsun_unit = h_cross_rsun./h_cross_rsun_norm;            % Unit n-vector 

% Projection of the Sun-Pointing vector onto the orbit plane and the satellite position vector
rsun_in_plane = cross(h_cross_rsun, orb_normal,2);              % Projection of the Sun-Pointing angle onto the orbital plane
rsun_in_plane_norm = sqrt(sum(abs(rsun_in_plane).^2,2));        % Magnitude of the Sun-pointing vector (for each column);
rsun_in_plane_norm = repmat(rsun_in_plane_norm,1,3);
rsun_in_plane_unit = rsun_in_plane./rsun_in_plane_norm;         % Unit vector of the projection of the Sun-Pointing angle onto the orbital plane   

% Angle ALPHA
num_alpha = dot(tx_ECI_pos_unit', h_cross_rsun_unit, 2);        % Numerator for the arctan to compute Alpha
den_aplha = dot(tx_ECI_pos_unit', rsun_in_plane_unit, 2);       % Denominator for the arctan to compute Alpha
alpha = atan2(num_alpha,den_aplha);                             % Angle ALPHA [rad]      

% Quadrant check: make sure 0 < alpha < 2*pi
alpha_corr = find(alpha<0);
alpha(alpha_corr) = alpha(alpha_corr) + 2*pi; 

%Nominal yaw
betar = deg2rad(beta);
for idx=1:length(betar)
    if isnan(betar(idx))
        yaw_all(idx) = NaN;
    elseif betar(idx)>0
        yaw_all(idx) = (atan(tan(betar(idx))/sin(alpha(idx))))*r2d;
        if (alpha(idx)>pi && alpha(idx)<=2*pi)
            yaw_all(idx) = yaw_all(idx)+180;
        end
    elseif betar(idx)<=0
        yaw_all(idx) = mod((atan(tan(betar(idx))/sin(alpha(idx))))*r2d,360);
        if (alpha(idx)> pi)
            yaw_all(idx) = yaw_all(idx)+180;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 'Noon/Midnight' and 'Shadow' Yaw Exclusion condition %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the vector from the satellite to the sun
% Note: need to reverse the satellite position vector so it points from satellite to earth
r_sun_sat = rsun_in_plane + (-1)*pos_vect';                      % Sat-Sun vector 
r_sun_sat_norm = sqrt(sum(abs(r_sun_sat).^2,2));                 % Magnitude of the Sun-Sat vector (for each column)
r_sun_sat_norm = repmat(r_sun_sat_norm,1,3);
r_sun_sat_unit = r_sun_sat./r_sun_sat_norm;                      % Unit vector of the Sat-Sun vector       

% Compute the Earth-Satellite-Sun angle (in radians)
dot_E = dot((-1)*tx_ECI_pos_unit', r_sun_sat_unit,2);
E = acos(dot_E);                                                 % Earth-Satellite-Sun angle [rad]

% Compute the actual yaw angle, B (in degrees), induced by the yaw bias, b, inserted in the
% satellite ACS
b = 0.5;                                                         % Yaw bias in degrees, set as of Nov. 1995
B = asin(0.0175*b*d2r./sin(E))*r2d;                              % Yaw Angle[deg]

% Compute nominal yaw rate
R = 0.135;     % maximal yaw rate, deg/sec,
% ftp://sideshow.jpl.nasa.gov/pub/GPS_yaw_attitude/nominal_yaw_rates
RR = 0.0017;   % maximal yaw rate rate, deg/sec^2 (Bar-Sever)

% "mu" is defined by Bar-Sever to be the orbit angle, measured from orbit
% midnight in the direction of motion. Thus, it is alpha + 180 deg. This
% angle is assumed to vary little over time, and its rate can be replaced
% by a constant.
mu_dot = 0.0083;        % deg/sec
b = 0.5;                % degrees, set as of Nov. 1995

% Compute B_dot (deg/sec), the rate of change of B
B_dot = -RR*b*cos(E).*cos(beta*d2r).*sin(alpha+pi)*mu_dot./(cos(B*d2r).*(sin(E)).^3);
% Compute the nominal yaw rate (deg/sec)
nom_yaw_rate = (mu_dot*tan(beta*d2r).*cos(alpha+pi)./(sin(alpha+pi).^2+tan(beta*d2r).^2)) + B_dot;

% Determine if satellite is in the "noon/midnight maneuver regime"; this occurs in
% the vicinity of orbit noon when the nominal yaw rate would be higher than
% the maximal yaw rate allowed by the satellite.
dt = 1;
ISNOON = zeros(length(time),1);
ISSHADOW = zeros(length(time),1);

for time_tag = 1:1:length(time)
    if time_tag == 1
        ISNOON(time_tag) = 0;
        ISSHADOW(time_tag) = 0;
        yaw(time_tag) = yaw_all(time_tag);
        yaw_rate(time_tag) = nom_yaw_rate(time_tag);
    end 
    
    if ((beta(time_tag) > -5) && (beta(time_tag) < 5) && (time_tag ~= 1))
        if ISNOON(time_tag) == 0
            % satellite is NOT in the 'noon/midnight' maneuver regime
            yaw(time_tag) = yaw_all(time_tag);
            yaw_rate(time_tag) = nom_yaw_rate(time_tag);
            
            if abs(nom_yaw_rate(time_tag)) > R
                % satellite is beginning to enter the 'noon/midnight' maneuver regime
                yaw(time_tag) = NaN;                  
                yaw_rate(time_tag) = R * sign(yaw_all(time_tag)-yaw_all(time_tag-1));
                ISNOON(time_tag) = 1;
            end
        else
            % satellite is already in 'noon/midnight' maneuver regime
            yaw(time_tag) = yaw(time_tag-1) + yaw_rate(time_tag-1) * dt;
            yaw_rate(time_tag) = R * sign(yaw_all(time_tag)-yaw_all(time_tag-1));
            
            if yaw(time_tag) > yaw_all(time_tag) && (yaw_all(time_tag-1)-yaw_all(time_tag))<=0
                % satellite is on the "upslope" and has overshot the nominal yaw angle;
                % satellite returns to nominal regime
                ISNOON(time_tag) = 0;
            elseif yaw(time_tag) < yaw_all(time_tag) && (yaw_all(time_tag-1)-yaw_all(time_tag))>=0
                % satellite is on the "downslope" and has undershot the nominal yaw angle;
                % satellite returns to nominal regime
                ISNOON(time_tag) = 0;
            end
        end   % end ISNOON loop
        
    else
        % beta not in 'noon/midnight' range; satellite is in nominal regime
        yaw(time_tag) = yaw_all(time_tag);
        yaw_rate(time_tag) = nom_yaw_rate(time_tag);
    end  % end beta loop


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 'Shadow' Yaw Exclusion condition %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnan(yaw(time_tag))
    % If SV is in 'shadow' regime, the Yaw get excluded
    shadow_flag = compute_shadow(r_sun_vector(time_tag,:), pos_vect(:,time_tag));
    if  shadow_flag == 1
%           yaw(time_tag) = NaN;    % Remove Yaw from pattern reconstruction
        ISSHADOW(time_tag) = 1;
    end
end % end if isnan


end % end for time

end % end function
