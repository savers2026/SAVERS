%--------------------------------------------------------------------------
% This Function calculates velocity of Point P(X';Y';Z') in Local ref.Sys
% taking into account the Earth Rotation in ECEF coord. and LLH ref.System
%--------------------------------------------------------------------------
% Author: Feb 14, 2010  by R.Giusto, r.giusto@gmail.com
% Update: Sep 2019 for the case of spaceborne GNSSR by Laura Dente (Tor Vergata University)
% New version: Feb 2022 rewritten ex novo by Davide Comite (Sapienza University)

% Purpose: This function calculates velocity of P(X'Y'Z') in local coord.
%          computing its ECEF and LLH position for taking into account
%          the Earth Rotation motion
%--------------------------------------------------------------------------

function [vLOCAL] = earth_rot2local2_GNSSR_mat(phi_ENU2local, pLocal, OenuINecef)

         Nx    = size(pLocal(:,:,1),2);
         Ny    = size(pLocal(:,:,1),1);

          pECEF = local2ecef_dc_mat(phi_ENU2local, pLocal, OenuINecef);
                     
         [lat, long, h] = xyz2llh_mat(pECEF);

         %-------------------------------------------------------------
         % velocity components [ VX, VY, VZ ] in ECEF reference system
         %-------------------------------------------------------------
         
         RX_radius_rot = sqrt(pECEF(:,:,1).^2 + pECEF(:,:,2).^2 + pECEF(:,:,3).^2);              % units are: [m]
         V_rot         = 2*pi./(24*3600).*RX_radius_rot; % units are: [m/s]
         clat          = cos(lat);
         slon          = sin(long);
         clon          = cos(long);
%         V_rot =0; %AMIR TEST
         vX            = -V_rot.*clat.*slon;      % units are: [m/s] in ECEF coord.system
         vY            =  V_rot.*clat.*clon;      % units are: [m/s] in ECEF coord.system

         %--------------------------------------------------------------------------------
         v_vectECEF(:,:,1) = vX;                  % speed [m/s] in ECEF coord. system 
         v_vectECEF(:,:,2) = vY;                  % speed [m/s] in ECEF coord. system 
         v_vectECEF(:,:,3) = zeros(Ny,Nx);        % speed [m/s] in ECEF coord. system 
         %--------------------------------------------------------------------------------
 
         v_vectECEF_IN(:,:,1) = v_vectECEF(:,:,1) + OenuINecef(1);
         v_vectECEF_IN(:,:,2) = v_vectECEF(:,:,2) + OenuINecef(2);
         v_vectECEF_IN(:,:,3) = v_vectECEF(:,:,3) + OenuINecef(3);
         
         %vLOCAL = ecef2local_dc_mat(phi_ENU2local, v_vectECEF + OenuINecef, OenuINecef);  %%[m/s]
         vLOCAL = ecef2local_dc_mat(phi_ENU2local, v_vectECEF_IN, OenuINecef);  %%[m/s]
                                                
end
