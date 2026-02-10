%% ************************************************************************
% MODULE NAME:      gen_noise_sys.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      D. Comite, N. Pierdicca
% @copyright:   Sapienza University of Rome
% Original version:      Oct 2021 by D. Comite, N. Pierdicca
% Main updates:          Oct 2022 by D. Comite
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% 
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [NoisySignalPowerAnt_LR, NoisySignalPowerAnt_RR, NLbb_LHCP, NLbb_RHCP, G_sys_lin_LHCP, G_sys_lin_RHCP, G_sys_lin_LHCP_nominal , G_sys_lin_RHCP_nominal] = gen_noise_sys(NoiseFreeSignalPowerAnt_LR, NoiseFreeSignalPowerAnt_RR, geoSYSp,Gain_LHCP_onsurf_lin,Gain_RHCP_onsurf_lin)


%% Set random seed for all random generator during process for Noise (for MM)
rng(geoSYSp.randseed)

Tcohe = geoSYSp.cohIntTime_sec;
Tinco = geoSYSp.incoIntTime_sec;
CRF   = geoSYSp.CodeRepetitionFreq;
Kboltzman = 1.380649e-23; % in J/K
Bandwidth = 1/Tcohe;

%% implementation losses
LN_lin_LHCP = geoSYSp.LN_lin_LHCP;
LN_lin_RHCP = geoSYSp.LN_lin_RHCP;
LS_lin_LHCP = geoSYSp.LS_lin_LHCP;
LS_lin_RHCP = geoSYSp.LS_lin_RHCP;

%% noise figure Front End
F0_LHCP   = geoSYSp.NoiseFigure_lin_LHCP;
F0_RHCP   = geoSYSp.NoiseFigure_lin_RHCP;

DeltaFdet = geoSYSp.DeltaF_lin;

%% Noise Figure
%gNF = geoSYSp.gNF;
mNF_LHCP  = geoSYSp.mNF_LHCP;
mNF_RHCP  = geoSYSp.mNF_RHCP;

%% coefficients Gain
%gG = geoSYSp.gG;
mG_LHCP  = geoSYSp.mG_LHCP;
mG_RHCP  = geoSYSp.mG_RHCP;

%% physical temperature of the receiver
TRx_K      = geoSYSp.TRx_K;
DeltaTRx_K = geoSYSp.DeltaTRx_K;

%% physical temperature with fluctuations in K
TRx_flu_K = TRx_K + DeltaTRx_K;

%% system gain
G_sys0_lin_LHCP = geoSYSp.G_sys0_lin_LHCP;
G_sys0_lin_RHCP = geoSYSp.G_sys0_lin_RHCP;
% 
% G_sys0_lin_LHCP = db2pow(33);
% G_sys0_lin_RHCP = db2pow(33);

DeltaG_sys_lin  = geoSYSp.DeltaG_sys_lin;

%% Total Gain with fluctuations
G_sys_lin_LHCP     = G_sys0_lin_LHCP.*10.^(mG_LHCP.*(TRx_flu_K - 273.15)./10);
G_sys_flu_lin_LHCP = G_sys_lin_LHCP.*DeltaG_sys_lin;
G_sys_lin_RHCP     = G_sys0_lin_RHCP.*10.^(mG_RHCP.*(TRx_flu_K - 273.15)./10);
G_sys_flu_lin_RHCP = G_sys_lin_RHCP.*DeltaG_sys_lin;
G_sys_lin_LHCP_nominal = LS_lin_LHCP.*G_sys0_lin_LHCP.*10.^(mG_LHCP.*(TRx_K - 273.15)./10);
G_sys_lin_RHCP_nominal = LS_lin_RHCP.*G_sys0_lin_RHCP.*10.^(mG_RHCP.*(TRx_K - 273.15)./10);

%G_sys_flu_dB = 10*log10(G_sys_lin.*DeltaG_sys_lin);

%% Noise Figure with fluctuations
F_LHCP         = F0_LHCP.*10.^(mNF_LHCP.*(TRx_flu_K - 273.15)./10);
F_flu_lin_LHCP = F_LHCP.*DeltaFdet;
F_RHCP         = F0_RHCP.*10.^(mNF_RHCP.*(TRx_flu_K - 273.15)./10);
F_flu_lin_RHCP = F_RHCP.*DeltaFdet;

if F_flu_lin_LHCP < 1 
    F_flu_lin_LHCP = 1; 
    disp(['******* Warning ******']), 
    disp(['The variation of the noise figure is bigger than its nominal value, F has been set to 0 dB'])
end
if F_flu_lin_RHCP < 1 
    F_flu_lin_RHCP = 1; 
    disp(['******* Warning ******']), 
    disp(['The variation of the noise figure is bigger than its nominal value, F has been set to 0 dB'])
end

%F_flu_dB  = 10*log10(F.*DeltaFdet);

%% Antenna temperature and receiver temperature
TA_K   = geoSYSp.TA_down_K;
Tref_K = geoSYSp.ReceiverTemperature_K;
% Nindependentsamples = Tinco*CRF;
Nindependentsamples = floor(Tinco/Tcohe); %AMIR TEST
%%

PN_LHCP = Kboltzman*(TA_K + (F_flu_lin_LHCP - 1)*Tref_K)*Bandwidth;
PN_RHCP = Kboltzman*(TA_K + (F_flu_lin_RHCP - 1)*Tref_K)*Bandwidth;

%%%DC: Noise Boost
PN_LHCP = PN_LHCP/10.^(geoSYSp.sigPowerBooster_dB/10);
PN_RHCP = PN_RHCP/10.^(geoSYSp.sigPowerBooster_dB/10);

PN_star_LHCP = LN_lin_LHCP.*G_sys_flu_lin_LHCP.*PN_LHCP;
PN_star_RHCP = LN_lin_RHCP.*G_sys_flu_lin_RHCP.*PN_RHCP;

% if there is LNA in the receiver RF-FRONT

Pr_star_LR = LS_lin_LHCP.*G_sys_flu_lin_LHCP.*NoiseFreeSignalPowerAnt_LR;
Pr_star_RR = LS_lin_RHCP.*G_sys_flu_lin_RHCP.*NoiseFreeSignalPowerAnt_RR;

% If only LNA is for active antenna  - AMIR
% Pr_star_LR = LS_lin_LHCP.*NoiseFreeSignalPowerAnt_LR;
% Pr_star_RR = LS_lin_RHCP.*NoiseFreeSignalPowerAnt_RR;

% rng(20);

nT_star_LR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);   
nT_star_RR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);  

nR_star_LR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);
nR_star_RR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);

NoiseSamples_LR = PN_star_LHCP.*(1 + nT_star_LR);
NoiseSamples_RR = PN_star_RHCP.*(1 + nT_star_RR);

%%%It is unlikely, but we can produce a negative noise sample the way
%%%it is generated above
NoiseSamples_LR(NoiseSamples_LR <0) = 0;
NoiseSamples_RR(NoiseSamples_RR <0) = 0;

speckleNoiseFactor_LR  = 1 + (1 - geoSYSp.alfa)*nR_star_LR;
speckleNoiseFactor_RR  = 1 + (1 - geoSYSp.alfa)*nR_star_RR;

SpeckleSamples_LR      = Pr_star_LR.*speckleNoiseFactor_LR;
SpeckleSamples_RR      = Pr_star_RR.*speckleNoiseFactor_RR;

% NoisySignalPowerAnt_LR = Pr_star_LR + NoiseSamples_LR;


NoisySignalPowerAnt_LR = SpeckleSamples_LR + NoiseSamples_LR;
NoisySignalPowerAnt_RR = SpeckleSamples_RR + NoiseSamples_RR;



%% --------------- Black Body DDM --------------------------------------

%load temperature with fluctuations
TL_K         = geoSYSp.TL_K + geoSYSp.DeltaTL_K;

%physical temperature  with fluctuations
DeltaTRxbb_K = geoSYSp.DeltaTRxbb_K;
TRxbb_flu    = TRx_K + DeltaTRxbb_K;

%noise Figure with fluctuations
DeltaFbb_lin      = geoSYSp.DeltaFbb_lin;
Fbb_lin_LHCP      = F0_LHCP.*10.^(mNF_LHCP.*(TRxbb_flu - 273.15)./10);
Fbb_flu_lin_LHCP  = Fbb_lin_LHCP.*DeltaFbb_lin;
Fbb_lin_RHCP      = F0_RHCP.*10.^(mNF_RHCP.*(TRxbb_flu - 273.15)./10);
Fbb_flu_lin_RHCP  = Fbb_lin_RHCP.*DeltaFbb_lin;

%gain fluctuation
DeltaG_sysbb_lin  = geoSYSp.DeltaG_sysbb_lin;
G_sys_flubb_LHCP  = G_sys_lin_LHCP.*DeltaG_sysbb_lin;
G_sys_flubb_RHCP  = G_sys_lin_RHCP.*DeltaG_sysbb_lin;

%PL     = Kboltzman*(TL + (10^(Fbb_flu/10) - 1)*Tref)*Bandwidth;
PL_LHCP    = Kboltzman*(TL_K + (Fbb_flu_lin_LHCP - 1)*Tref_K)*Bandwidth;
PL_RHCP    = Kboltzman*(TL_K + (Fbb_flu_lin_RHCP - 1)*Tref_K)*Bandwidth;

nL_star_LR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);
nL_star_RR = randn(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay)./sqrt(Nindependentsamples);

PL_star_LHCP = PL_LHCP.*G_sys_flubb_LHCP.*LN_lin_LHCP;
NLbb_LHCP    = PL_star_LHCP.*(1 + nL_star_LR);
PL_star_RHCP = PL_RHCP.*G_sys_flubb_RHCP.*LN_lin_RHCP;
NLbb_RHCP    = PL_star_RHCP.*(1 + nL_star_RR);

end