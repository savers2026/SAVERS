%--------------------------------------------------------------------------
% This Function calculates DOPPLER shift for Scattering Points(X';Y';Z') 
% belonging to the Area of Interest in the Local Reference System
%--------------------------------------------------------------------------
% Author: Feb 14, 2010  by R.Giusto, r.giusto@gmail.com
% Update: Sep 2019 for the case of spaceborne GNSSR by Laura Dente (Tor
%           Vergata University)
%
%          ###################################  
% Purpose: This function calculates doppler of P(X'Y'Z') in local coord.
%          taking into account all possible contributions
% Input :  xp, yp, zp    [2] X'Y'Z' position vector in [m] Local Coord. System
%          Tx        [3] position of TxGPS in    [m]
%          VTX_local [3] velocity of TxGPS in  [m/s] in Local System
%          Rx        [3] position of Rx    in    [m]
%          VRx_local [3] velocity of Rx    in  [m/s] in Local System
% Output:  DDR       [1] Doppler Shift in [Hz]
%
% Subroutine: earth_rot2local2_GNSSR
%--------------------------------------------------------------------------
% inner product :  u (row vector)   v (column vector)
%                  u = [3; 1; 4];
%                                   v = [2 0 -1];
%                  x = v*u =  2
%--------------------------------------------------------------------------

function [DDR,VRx_local] = refl_doppler2_GNSSR(phi_ENU2local, xp, yp, zp, Tx, VTX_local, Rx, VRx_local, OenuINecef, WaveL,TRAp) 

         %% AMIR: Adding Earth rotation instead of applying it on the Flight track simulation - GREAT
         if TRAp.Rxlocal(3) < 15e3
          vRxECEF_EARTH_ROT = TRAp.VRxECEF + (cross([0,0,7.2921158553*1.e-5],TRAp.RxECEF)).';
          vRxENU_EARTH_ROT = ecef2enu_dc(vRxECEF_EARTH_ROT,TRAp.OENUinECEF, 'v');   
          VRx_local  = rotation_clockwise(vRxENU_EARTH_ROT,   TRAp.phi_ENU2local, 'rad');
%           vHydroG_ECEF  = vHydroG_ECEF + (cross([0,0,7.2921158553*1.e-5],pHydroG_ECEF)).';
         end
         %% 


        P_in  = [xp ; yp ; zp];          % units are [m]
        vect  = P_in       ;             % row vector            
        V_loc = earth_rot2local2_GNSSR(phi_ENU2local, P_in, OenuINecef);  % [m/s]

        %versor (Rx)-to-(P_in) is :
        %----------------------------
        
        u    = vect-Rx;
        v    = transpose(u);
        modu = sqrt(v*u);
        uRP  = u / modu;
        vRP  = transpose(uRP);

        %
        %versor (Tx)-to-(P_in) is :
        %----------------------------

        u    = Tx-vect;
        v    = transpose(u);
        modu = sqrt(v*u);
        uPT  = u / modu;
        vPT  = transpose(uPT);
 
       %DOPPLER COMPUTATION :
       % -->  [Hertz]
       
        DDR = vPT*(V_loc - VTX_local)/WaveL + vRP*(VRx_local - V_loc)/WaveL;

        %--------------------------------------------------------------------------
end