%--------------------------------------------------------------------------
% This Function calculates DOPPLER shift for Scattering Points(X';Y';Z') 
% belonging to the Area of Interest in the Local Reference System
%--------------------------------------------------------------------------
% Author: Feb 14, 2010  by R.Giusto, r.giusto@gmail.com
% Update: Sep 2019 for the case of spaceborne GNSSR by Laura Dente (Tor Vergata University)
% New version: Feb 2022 rewritten ex novo by Davide Comite (Sapienza University)
%
%          ###################################  
% Purpose: This function calculates doppler of the point on surface in local coord.
%          taking into account all possible contributions

% Subroutine: earth_rot2local2_GNSSR

function DDR = refl_doppler2_GNSSR_mat(pSGR, TRAp, WaveL) 




        
         phi_ENU2local = TRAp.phi_ENU2local;  
         Tx            = TRAp.Txlocal;
         VTx_local     = TRAp.VTxlocal;
         Rx            = TRAp.Rxlocal; 
         VRx_local     = TRAp.VRxlocal;
         OenuINecef    = TRAp.OENUinECEF;

         Nx    = size(pSGR(:,:,1),2);
         Ny    = size(pSGR(:,:,1),1);        

         %% AMIR: Adding Earth rotation instead of applying it on the Flight track simulation - GREAT
         if TRAp.Rxlocal(3) < 15e3
          vRxECEF_EARTH_ROT = TRAp.VRxECEF + (cross([0,0,7.2921158553*1.e-5],TRAp.RxECEF)).';
          vRxENU_EARTH_ROT = ecef2enu_dc(vRxECEF_EARTH_ROT,TRAp.OENUinECEF, 'v');   
          VRx_local  = rotation_clockwise(vRxENU_EARTH_ROT,   TRAp.phi_ENU2local, 'rad');
%           vHydroG_ECEF  = vHydroG_ECEF + (cross([0,0,7.2921158553*1.e-5],pHydroG_ECEF)).';
         end
         %% 

         V_loc = earth_rot2local2_GNSSR_mat(phi_ENU2local, pSGR, OenuINecef);  % [m/s]

        
         %Unit vector Rx to P_in         
         ux    = pSGR(:,:,1) - Rx(1);
         uy    = pSGR(:,:,2) - Rx(2);
         uz    = pSGR(:,:,3) - Rx(3);

         modu  = sqrt(ux.^2 + uy.^2 + uz.^2);
                  
         uxRP  = ux./modu;
         uyRP  = uy./modu;
         uzRP  = uz./modu;

         %Unit vector Tx to P_in
         vx    = Tx(1) - pSGR(:,:,1);
         vy    = Tx(2) - pSGR(:,:,2);
         vz    = Tx(3) - pSGR(:,:,3);         

         modu  = sqrt(vx.^2 + vy.^2 + vz.^2);
         
         vxPT  = vx./modu;
         vyPT  = vy./modu;
         vzPT  = vz./modu;

         %VTX_x_local = repmat(VTX_local(1), Ny, Nx);
         %VTX_y_local = repmat(VTX_local(2), Ny, Nx);
         %VTX_z_local = repmat(VTX_local(3), Ny, Nx);

         %DOPPLER COMPUTATION  % -->  [Hertz]       
         %DDR = vPT*(V_loc - VTx_local)/WaveL + vRP*(VRx_local - V_loc)/WaveL;

         DDR = (vxPT.*(V_loc(:,:,1) - VTx_local(1)) + vyPT.*(V_loc(:,:,2) - VTx_local(2)) + vzPT.*(V_loc(:,:,3) - VTx_local(3)))./WaveL +...
               (uxRP.*(VRx_local(1) - V_loc(:,:,1)) + uyRP.*(VRx_local(2) - V_loc(:,:,2)) + uzRP.*(VRx_local(3) - V_loc(:,:,3)))./WaveL;

         %--------------------------------------------------------------------------
end