%--------------------------------------------------------------------------
% This Function calculates DOPPLER shift for Scattering Points(X';Y';Z') 
% belonging to the Area of Interest in the Local Reference System
%--------------------------------------------------------------------------
function DDR = refl_doppler2(xp, yp, Tx, VTX_local, Rx, VRx_local) 

% last revision: Feb 14, 2010  by R.Giusto
% r.giusto@gmail.com
% rgiusto@libero.it
%
%             Function DDR = refl_doppler( )
%          ###################################  
% Purpose: This function calculates doppler of P(X'Y'Z') in local coord.
%          taking into account all possible contributions
% Input :  xp, yp    [2] X'Y' position vector in [m] Local Coord. System
%          Tx        [3] position of TxGPS in    [m]
%          VTX_local [3] velocity of TxGPS in  [m/s] in Local System
%          Rx        [3] position of Rx    in    [m]
%          VRx_local [3] velocity of Rx    in  [m/s] in Local System
% Output:  DDR       [1] Doppler Shift in [Hz]
%--------------------------------------------------------------------------
% inner product :  u (row vector)   v (column vector)
%                  u = [3; 1; 4];
%                                   v = [2 0 -1];
%                  x = v*u =  2
%--------------------------------------------------------------------------
WaveL= 0.299792458/1.5754;                   % [m]wavelenght of GPS signal
%speedlight = 2.99792458e2;                  % [m/usec]
        P_in = [xp ; yp ; 0];            % units are [m]
        vect = P_in       ;              % row vector
        V_loc = [0;0;0]   ;             
        V_loc = earth_rot2local2(P_in);  % [m/s]
       
% versor (Rx)-to-(P_in) is :
%----------------------------
%
        u=Rx-vect;
        v=transpose(u);
        modu=sqrt(v*u);
        uRP = u / modu;
        vRP = transpose(uRP);
%
% versor (Tx)-to-(P_in) is :
%----------------------------
%
        u=Tx-vect;
        v=transpose(u);
        modu=sqrt(v*u);
        uTP = u / modu;
        vTP = transpose(uTP);
%
%       DOPPLER COMPUTATION :
%---------------------------->  [Hertz]
%
        DDR = vRP*(VRx_local-V_loc)/WaveL + vTP*(VTX_local-V_loc)/WaveL ;
       
%--------------------------------------------------------------------------
end