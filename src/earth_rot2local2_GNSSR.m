%--------------------------------------------------------------------------
% This Function calculates velocity of Point P(X';Y';Z') in Local ref.Sys
% taking into account the Earth Rotation in ECEF coord. and LLH ref.System
%--------------------------------------------------------------------------
% Author: Feb 14, 2010  by R.Giusto, r.giusto@gmail.com
% Update: Sep 2019 for the case of spaceborne GNSSR by Laura Dente (Tor Vergata University)
% Update: Nov 2021 for HSAVES - Davide Comite

% Purpose: This function calculates velocity of P(X'Y'Z') in local coord.
%          computing its ECEF and LLH position for taking into account
%          the Earth Rotation motion
% Input :  P_in    [3] X'Y'Z' position vector in Local Coord. System
% Output:  V_local [3] VX'VY'VZ' velocity vector in the Local System
% Subroutine: xyz2llh, ecef2local2_GNSSR
%--------------------------------------------------------------------------

function [vLOCAL] = earth_rot2local2_GNSSR(phi_ENU2local, pLocal, OenuINecef)

          pECEF = local2ecef_dc(phi_ENU2local, pLocal, OenuINecef);
           
          X = pECEF(1);
          Y = pECEF(2);
          Z = pECEF(3);
          
         [lat, long, h] = xyz2llh(X,Y,Z);

         %--------------------------------------------------------------------------
         % velocity components [ VX, VY, VZ ] in ECEF reference system
         %-------------------------------------------------------------------
         
         RX_radius_rot = sqrt(X^2 + Y^2 + Z^2);              % units are: [m]
         V_rot         = 2*pi/(24*3600)*RX_radius_rot; % units are: [m/s]
         clat          = cos(lat);
         slon          = sin(long);
         clon          = cos(long);
%         V_rot =0; %AMIR TEST
         vX            = -V_rot*clat*slon;             % units are: [m/s] in ECEF coord.system
         vY            =  V_rot*clat*clon;             % units are: [m/s] in ECEF coord.system

         %--------------------------------------------------------------------------
         v_vectECEF    = [vX; vY; 0];                  % speed [m/s] in ECEF coord. system 
         %--------------------------------------------------------------------------
 
         
         vLOCAL = ecef2local_dc(phi_ENU2local, v_vectECEF + OenuINecef, OenuINecef);  %%[m/s]
                                                
end
