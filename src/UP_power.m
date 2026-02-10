%% ************************************************************************
% MODULE NAME:      UP_power_routine.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      D. Comite and L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      2020 by D. Comite and L. Guerriero
% Main updates:          May 2022 by D. Comite; Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% 
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [LR_UP_Power, RR_UP_Power, LR_UP_Power_Flu, RR_UP_Power_Flu, PN_star, NLbb_UP, UP_pattern,Transm_pattern] = UP_power(TRAp, geoSYSp, indexTxFrequency)

GPS_theta         = TRAp.theta_Txpattern;
GPS_phi           = TRAp.phi_Txpattern; 
GPSPattern_linear = TRAp.GainTx_linear(:,:,indexTxFrequency);

incAngle_antennaPattern = TRAp.theta_RxUPpattern;
azimuthAngle_antennaPattern  = TRAp.phi_RxUPpattern; 
antennaPattern_linear   = TRAp.GainRxUP(:,:,indexTxFrequency);

%% Distances

RTx       = (TRAp.TxECEF - TRAp.RxECEF);

% distance in [m] GPS_Tx - Rx ECEF
RTx_mod   = norm(RTx);

factor_UP = (TRAp.waveL(indexTxFrequency)./(4*pi*RTx_mod)).^2; %Original
% factor_UP2 =geoSYSp.cohIntTime_sec^2*(TRAp.waveL./(4*pi*RTx_mod)).^2;


%% Transform UP Rx from ECEF in Antenna Frame

if geoSYSp.Rx_Attitude_Variation_RPY(1) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(2) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(3) == 1
RxAttitude = [TRAp.Rx_roll,TRAp.Rx_pitch,0];
else
RxAttitude =  TRAp.Rx_nomAtt_angles;
end

RxAttitude = -1 .* RxAttitude;

[Rxmat_ECEF2OF, Rxmat_OF2BF, Rxmat_BFtoAF] = matrixFrameGen(TRAp.RxECEF, TRAp.VRxECEF, RxAttitude, TRAp.RxAttitude_angles, TRAp.RxEuler_angles);

Rxmat_BFtoAF(:,2) = -Rxmat_BFtoAF(2,:);
%UP antenna is looking zenith
Rxmat_BFtoAF(:,3) = -Rxmat_BFtoAF(3,:);

RTx_OF = Rxmat_ECEF2OF*RTx;
RTx_BF = Rxmat_OF2BF*RTx_OF;
RTx_AF = Rxmat_BFtoAF*RTx_BF;

%piercing angles
thetaUP_pier  = acosd(dot(RTx_AF./norm(RTx_AF), [0,0,1]));
phiUP_pier    = atan2d(RTx_AF(2), RTx_AF(1));

[~,ind_phiUP] = min(abs(phiUP_pier   - azimuthAngle_antennaPattern));
[~,ind_thUP]  = min(abs(thetaUP_pier - incAngle_antennaPattern));

%% Transform Tx in Antenna Frame

%correcting yaw of the TX
TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles(3) +  TRAp.Yaw_TX;

[Txmat_ECEF2OF, Txmat_OF2BF, Txmat_BFtoAF] = matrixFrameGen(TRAp.TxECEF,  TRAp.VTxECEF, TRAp.Tx_nomAtt_angles, TRAp.TxAttitude_angles, TRAp.TxEuler_angles);

TRx_K = (TRAp.RxECEF - TRAp.TxECEF);

TRx_OF = Txmat_ECEF2OF*TRx_K;
TRx_BF = Txmat_OF2BF*TRx_OF;
TRx_AF = Txmat_BFtoAF*TRx_BF;

%piercing angles
thetaGPS_pier = acosd(dot(TRx_AF./norm(TRx_AF), [0,0,1]));
phiGPS_pier   = atan2d(TRx_AF(2), TRx_AF(1));

%DC: comment added on oct 2023
% [~,ind_phiTransm] = min(abs(abs(phiGPS_pier) - GPS_phi)); 
[~,ind_phiTransm] = min(abs(phiGPS_pier - GPS_phi));
[~,ind_thTransm]  = min(abs(thetaGPS_pier    - GPS_theta));

UP_pattern     = antennaPattern_linear(ind_phiUP, ind_thUP);
Transm_pattern = GPSPattern_linear(ind_phiTransm, ind_thTransm);

[pt, prR, prL] = polarization_unit_vectors(TRAp.flagTxMismatch, TRAp.flagRxUPmismatch, thetaGPS_pier, thetaUP_pier, 'deg');

pt_hat  = pt./norm(pt);
prR_hat = prR./norm(prR);
prL_hat = prL./norm(prL);

PLF_up_RR = dot(prR_hat,pt_hat)*conj(dot(prR_hat,pt_hat));
PLF_up_LR = dot(prL_hat,pt_hat)*conj(dot(prL_hat,pt_hat));

LR_UP_Power = factor_UP.*PLF_up_LR.*UP_pattern.*Transm_pattern*db2pow(TRAp.TxSignalPower(indexTxFrequency));
RR_UP_Power = factor_UP.*PLF_up_RR.*UP_pattern.*Transm_pattern.*db2pow(TRAp.TxSignalPower(indexTxFrequency));


%% fluctuation on the direct signal power: defintion of variables

% defining variables
LN_lin = geoSYSp.LN_lin_RHCP_zenith;
LS_lin = geoSYSp.LS_lin_RHCP_zenith;

% coefficients Noise Figure
%gNF = geoSYSp.gNF;
mNF = geoSYSp.mNF_RHCP_zenith;

% coefficients Gain
%gG = geoSYSp.gG;
mG = geoSYSp.mG_RHCP_zenith;

Tcohe = geoSYSp.cohIntTime_sec;
Tinco = geoSYSp.incoIntTime_sec;
CRF   = geoSYSp.CodeRepetitionFreq;

% Nindependentsamples = Tinco*CRF;
Nindependentsamples = floor(Tinco/Tcohe); %AMIR TEST


Kboltzman = 1.380649e-23; % in J/K
Bandwidth = 1/Tcohe;

% physical temperature of the receiver
TRx_K       = geoSYSp.TRx_K;
DeltaTRx_up = geoSYSp.DeltaTRx_up_K;

TA_K  = geoSYSp.TA_up_K;
Tref_K = geoSYSp.ReceiverTemperature_K;

F0_lin        = geoSYSp.NoiseFigure_lin_RHCP_zenith;
DeltaF_up_lin = geoSYSp.DeltaF_up_lin;

% [Ndelay, Ndoppler]  = size(NoiseFreeSignalPowerAnt_LR);

%% fluctuation direct signal power

% physical temperature  with fluctuations
TRx_flu_K = TRx_K + DeltaTRx_up;

% system gain
G_sys0_lin       = geoSYSp.G_sys0_lin_RHCP_zenith;
DeltaG_sys_uplin = geoSYSp.DeltaG_sys_up_lin;

% Total Gain with fluctuations
G_sys_lin     = G_sys0_lin.*10.^(mG.*(TRx_flu_K - 273)./10);
G_sys_flu_lin = G_sys_lin.*DeltaG_sys_uplin; 

Pd_star_LR = LS_lin.*G_sys_flu_lin.*LR_UP_Power;
Pd_star_RR = LS_lin.*G_sys_flu_lin.*RR_UP_Power;

% Noise Figure with fluctuations
F         = F0_lin.*10.^(mNF.*(TRx_flu_K - 273)./10);
F_flu_lin = F.*DeltaF_up_lin;
%F_flu_dB = 10*log10(F.*DeltaF_up);

PN              = Kboltzman*(TA_K + (F_flu_lin - 1)*Tref_K)*Bandwidth;

%%%DC: Noise Boost
PN = PN/10.^(geoSYSp.sigPowerBooster_dB/10);

PN_star         = LN_lin.*G_sys_flu_lin.*PN;
% rng(20);

% nT_star       = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);
%no need of DDM BB for direct signal
nT_star         = randn(1)./sqrt(Nindependentsamples);
NoiseSamples    = PN_star.*(1 + nT_star);

LR_UP_Power_Flu = Pd_star_LR + NoiseSamples;
RR_UP_Power_Flu = Pd_star_RR + NoiseSamples;

% %% === AMIRREZA TEST X
% 
% % constants (unchanged)
% BW_RF_Hz   = 2.5e6;
% BW_DEC_Hz = 32e3;
% Tcohe      = 0.001;  % 1e-3
% Tinco      = 1; % 1.0
% Ncoh       = round(BW_DEC_Hz*Tcohe);  % 32
% Nlooks     = Tinco/Tcohe;             % 1000
% 
% % Stage-0: RF‐Analog
% Sig_RF_W     = RR_UP_Power;
% Noise_RF_W   = Kboltzman*(TA_K + (F_flu_lin - 1)*Tref_K)*BW_RF_Hz;
% 
% % Stage-1: Pre-Corr
% Sig_pre_W    = Pd_star_RR;
% Noise_pre_W  = LN_lin.*G_sys_flu_lin.*Noise_RF_W;
% 
% % Stage-2: Post-Coherent
% Sig_postC_cnt   = Sig_pre_W   * Ncoh^2;
% Noise_postC_cnt = Noise_pre_W * sqrt(Ncoh);
% 
% % Stage-3: Post-Integrate
% Sig_postI_cnt   = Sig_postC_cnt;
% Noise_postI_cnt = Noise_postC_cnt/sqrt(Nlooks);
% Sum_Sig_Noise_int_amir = Sig_postI_cnt + Noise_postI_cnt;
% 
% amir_out_test.Sig_postI = Sig_postI_cnt;
% amir_out_test.Noise_postI_cnt = Noise_postI_cnt;
% amir_out_test.Sum_Sig_Noise_int_amir = Sum_Sig_Noise_int_amir;
% amir_out_test.RR_UP_Power_Flu_def = RR_UP_Power_Flu;
% amir_out_test.NoiseSamples_def = NoiseSamples;
% amir_out_test.RR_UP_Power_def = RR_UP_Power;
% 
% %  %% ==== End Amirreza Post correlator ====

%% --------------- Black Body DDM --------------------------------------

%load temperature with fluctuations
TL_K      = geoSYSp.TL_K + geoSYSp.DeltaTL_K;

%physical temperature  with fluctuations
DeltaTRxbb_up = geoSYSp.DeltaTRxbb_up;
TRxbb_flu     = TRx_K + DeltaTRxbb_up;

%noise Figure with fluctuations
DeltaFbb_up_lin = geoSYSp.DeltaFbb_up_lin;
Fbb_lin         = F0_lin.*10.^(mNF.*(TRxbb_flu - 273)./10);
Fbb_flu_lin     = Fbb_lin.*DeltaFbb_up_lin;
%Fbb_flu_dB     = 10*log10(Fbb.*DeltaFbb_up);

%gain fluctuation
DeltaG_sysbb_lin = geoSYSp.DeltaG_sysbb_up_lin;
G_sys_flubb_lin  = G_sys_lin.*DeltaG_sysbb_lin;

%PL     = Kboltzman*(TL + (10^(Fbb_flu/10) - 1)*Tref)*Bandwidth;
PL      = Kboltzman*(TL_K + (Fbb_flu_lin - 1)*Tref_K)*Bandwidth;

% need for a DDM for the UP BB
nL_star = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);
% no need of DDM BB for direct signal
% nL_star = randn(1)./sqrt(Nindependentsamples);
PL_star = PL.*G_sys_flubb_lin.*LN_lin;
NLbb_UP = PL_star.*(1 + nL_star);


end