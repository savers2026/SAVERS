%% ************************************************************************
% MODULE NAME:      TxGain_SP.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      D. Comite
% @copyright:   Sapienza University of Rome
% Original version:      2021 by D. Comite
% Main updates:          Oct 2022 by D. Comite; Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% New updates:           Apr 2023 by D. Comite (Yaw Steering)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module computes the Tx antenna gain at the specular point (after pearcing).
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [Yaw_TX, GPSgainAtSP,TxAttitude_angles_yaw] = TxGain_SP(TRAp, GNSS_Week, GNSS_SoW)


%% Yaw correction
% posECEF = (3xn) position vector of the satellite (m, ECI)
% velECEF = (3xn) velocity vector of the satellite (m/s, ECI)

% gpsWeek = week number
gpsWeek      = GNSS_Week;
% gpsSeconds = time in seconds
gpsSeconds   = GNSS_SoW;
% gpsUTCOffset = leap seconds
gpsUTCOffset = 0;
% Calculate Julian date
julianDate   = 2444244.5 + (gpsWeek * 7) + ((gpsSeconds - gpsUTCOffset) / 86400.0);

% Calculate sidereal time
thetaGMST = lstime(julianDate);

% ECEF to ECI

posECEF = TRAp.TxECEF;
velECEF = TRAp.VTxECEF;

[posECI, velECI]= ecef2eci(posECEF, thetaGMST, velECEF);

TXposECI = posECI;
TXvelECI = velECI;

% Calculate yaw

T_JC       = (julianDate - 2451545)/36525;
TX_dir_ECI = TXposECI./norm(TXposECI); % SV RX unit direction vector

% [Yaw_TX, Alpha_TX, Beta_TX] = get_beta_yaw_angles_old(TXposECI, TXvelECI, T_JC, gpsSeconds, TX_dir_ECI)
[Yaw_TX, Alpha_TX, Beta_TX] = get_beta_yaw_angles(TXposECI, TXvelECI, T_JC, gpsSeconds, TX_dir_ECI);

%correcting yaw of the TX
TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles(3) +  Yaw_TX;
TxAttitude_angles_yaw = TRAp.TxAttitude_angles(3);

%% Transform Tx in Antenna Frame

%remind OENUinECEF is the specular point in ECEF
TxSP = (TRAp.OENUinECEF - TRAp.TxECEF);

[Txmat_ECEF2OF, Txmat_OF2BF, Txmat_BFtoAF] = matrixFrameGen(TRAp.TxECEF, TRAp.VTxECEF, TRAp.Tx_nomAtt_angles, TRAp.TxAttitude_angles, TRAp.TxEuler_angles);

TxSP_OF = Txmat_ECEF2OF*TxSP;
TxSP_BF = Txmat_OF2BF*TxSP_OF;
TxSP_AF = Txmat_BFtoAF*TxSP_BF;

%piercing angles
thetaGPS_pier = acosd(dot(TxSP_AF./norm(TxSP_AF), [0,0,1]));
phiGPS_pier   = atan2d(TxSP_AF(2), TxSP_AF(1));

%[val,ind_phiTransm] = min(abs(abs(phiGPS_pier) - TRAp.phi_Txpattern));
[val,ind_phiTransm] = min(abs(phiGPS_pier      - TRAp.phi_Txpattern));
[val,ind_thTransm]  = min(abs(thetaGPS_pier    - TRAp.theta_Txpattern));

GPSgainAtSP(:) = TRAp.GainTx_linear(ind_phiTransm,ind_thTransm,:);


end