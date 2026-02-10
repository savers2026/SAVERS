%% ************************************************************************
% MODULE NAME:      BISTATIC_IHS_dc.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      D. Comite, L. Dente, N. Pierdicca, L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      2010 by Roberto Giusto (airborne)
% Main updates:          2018 by L. Dente (spaceborne), 
%                        2021 D. Comite (HydroGNSS)
%                        2022 D. Comite and L.Dente
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel, config. param. for Lnadir&Rnadir&zenith)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module integrates the scattering of all facets in the simulation area 
% computing the power DDM free of noise and with noise and calls the routines to compute 
% the coherent channel, the thermal and spekle noise, the receiver noise, 
% the zenith antenna power
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [max_sigPower_noNoise_LR, max_sigPower_noNoise_RR, max_sigPower_Noise_LR, max_sigPower_Noise_RR,...
          sigPower_noNoise_LR_lin, sigPower_noNoise_RR_lin, sigPower_Noise_LR_lin, sigPower_Noise_RR_lin, NLbb_LR, NLbb_RR,...
          LR_UP_Power, RR_UP_Power, LR_UP_Power_Flu, RR_UP_Power_Flu, PN_star, NLbb_UP, UP_pattern, TxGainAtRx, geoSYSp,...
          AcohLR, AdiffLR, rho_TRsp] = BISTATIC_IHS_dc(Ny_surf, Nx_surf, sigmaRR_cohe, sigmaRR_inco, sigmaLR_cohe, sigmaLR_inco,...         
                                                       DDMparam, geoSYSp, TRAp, nameOutputFile, tag, surfParam, instrumentConfParam, indexTxFrequency)
            
            
%% MAPPING in the Doppler/Delay plane ---------------------------
flagRxmismatch = TRAp.flagRxDOWNmismatch; % added by amir

geoSYSp.chipTime = [1.0 0.1]; %[L1/E1 L5/E5] chip time in microsec [usec]

% Freq       = (0:geoSYSp.ny_Dfreq - 1).*geoSYSp.fpas - geoSYSp.FRQmax;
Freq       = (0:geoSYSp.ny_Dfreq - 1).*geoSYSp.fpas - (geoSYSp.FRQmax/2 - geoSYSp.fpas/2);  % Added by amir for Rongowai
Freq_3D    = reshape(Freq,1,1,geoSYSp.ny_Dfreq);
Freq_array = repmat(Freq_3D,Ny_surf,Nx_surf,1);

argsinc = geoSYSp.cohIntTime_sec*(Freq_array - DDMparam.fdoppler(:,:,indexTxFrequency));    
sincF   = sinc_dc(argsinc);  

sigPower_noNoise_RR_cohe = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay);
sigPower_noNoise_RR_inco = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay);
sigPower_noNoise_LR_cohe = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay);
sigPower_noNoise_LR_inco = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay);


for lmap = 1 : geoSYSp.nx_Rdelay

    tau     = zeros(geoSYSp.nx_Rdelay,1);
        
    %delay
    tau(lmap) = (lmap-1) * geoSYSp.mpas(indexTxFrequency) - geoSYSp.delayOffset_microsec(indexTxFrequency);       
    % tau(lmap) = (lmap-1) * geoSYSp.mpas(indexTxFrequency) + DDMparam.absoluteDelayAtSP_microsec;       %Test AMIR

    ADR = abs(tau(lmap) - DDMparam.fdelay(:,:,indexTxFrequency)); % <<<<<<<< [usec]

    %%

    modulationTxSignal = zeros(Ny_surf, Nx_surf);  
        
    conditionParam = ADR./geoSYSp.chipTime(indexTxFrequency);

    if geoSYSp.indexTxSignal==2 || geoSYSp.indexTxSignal==3 && indexTxFrequency==1         
        i0p5    = conditionParam <= 0.5;
        i0p51p0 = conditionParam > 0.5 & conditionParam <= 1;             
        modulationTxSignal(i0p5)    = 1-3*ADR(i0p5)./geoSYSp.chipTime(indexTxFrequency);
        modulationTxSignal(i0p51p0) = ADR(i0p51p0)./geoSYSp.chipTime(indexTxFrequency) - 1;
    else
        i1p0 = conditionParam <= 1;
        modulationTxSignal(i1p0) = 1-ADR(i1p0)/geoSYSp.chipTime(indexTxFrequency);
    end

        addendum_first = (sincF.*DDMparam.PATH_factor.*modulationTxSignal).^2;      %1 msec coherent integration time

    if flagRxmismatch == 0
        
    addendum_ci_LR = addendum_first.*DDMparam.Gain_LHCP_onsurf_lin(:,:,indexTxFrequency).*1;
    addendum_ci_RR = addendum_first.*DDMparam.Gain_RHCP_onsurf_lin(:,:,indexTxFrequency).*1;

    SUMRR_cohe = addendum_ci_RR.*sigmaRR_cohe.*DDMparam.factor(:,:,indexTxFrequency);
    SUMLR_cohe = addendum_ci_LR.*sigmaLR_cohe.*DDMparam.factor(:,:,indexTxFrequency);
    SUMRR_inco = addendum_ci_RR.*sigmaRR_inco.*DDMparam.factor(:,:,indexTxFrequency);
    SUMLR_inco = addendum_ci_LR.*sigmaLR_inco.*DDMparam.factor(:,:,indexTxFrequency);
    
    end
    % added by Amir , In order to have MM based on the cross-pol antenna
    % gain and not by polarization unit vector
    if flagRxmismatch == 1

    addendum_ci_LR = addendum_first;
    addendum_ci_RR = addendum_first;

    SUMRR_cohe = addendum_ci_RR.*((sigmaRR_cohe.*DDMparam.Gain_RHCP_onsurf_lin)+(sigmaLR_cohe.*DDMparam.Gain_RL_onsurf_lin)).*DDMparam.factor(:,:,indexTxFrequency);
    SUMLR_cohe = addendum_ci_LR.*((sigmaLR_cohe.*DDMparam.Gain_LHCP_onsurf_lin)+(sigmaRR_cohe.*DDMparam.Gain_LR_onsurf_lin)).*DDMparam.factor(:,:,indexTxFrequency);
    SUMRR_inco = addendum_ci_RR.*((sigmaRR_inco.*DDMparam.Gain_RHCP_onsurf_lin)+(sigmaLR_inco.*DDMparam.Gain_RL_onsurf_lin)).*DDMparam.factor(:,:,indexTxFrequency);
    SUMLR_inco = addendum_ci_LR.*((sigmaLR_inco.*DDMparam.Gain_LHCP_onsurf_lin)+(sigmaRR_inco.*DDMparam.Gain_LR_onsurf_lin)).*DDMparam.factor(:,:,indexTxFrequency);

    end

    SUMRR_cohe = squeeze(sum(sum(SUMRR_cohe,'omitnan')));
    SUMLR_cohe = squeeze(sum(sum(SUMLR_cohe,'omitnan')));
    SUMRR_inco = squeeze(sum(sum(SUMRR_inco,'omitnan')));
    SUMLR_inco = squeeze(sum(sum(SUMLR_inco,'omitnan')));
        
    %signal power down-looking antenna, no noise, no calibration
    sigPower_noNoise_RR_cohe(:, lmap) = SUMRR_cohe.';
    sigPower_noNoise_RR_inco(:, lmap) = SUMRR_inco.';
    sigPower_noNoise_LR_cohe(:, lmap) = SUMLR_cohe.';
    sigPower_noNoise_LR_inco(:, lmap) = SUMLR_inco.';

end

%% Power on surface for Covariance

%%%without areas - wrong
% DDMparam.PonSurf = DDMparam.factor./DDMparam.regGridSim_cellareas.*DDMparam.Gain_LHCP_onsurf_lin.*(sigmaLR_cohe + sigmaLR_inco).*DDMparam.PATH_factor^2;

%%% DDMparam.factor = TRAp.waveL^2*db2pow(TRAp.TxSignalPower)*TRAp.TxGainAtSP*regGridSim_cellareas./(4*pi)^3;
DDMparam.PonSurf = DDMparam.factor(:,:,indexTxFrequency)./1.*DDMparam.Gain_LHCP_onsurf_lin(:,:,indexTxFrequency).*...
                       (sigmaLR_cohe + sigmaLR_inco).*DDMparam.PATH_factor^2;

%% total LR and RR signal power, no noise, no calibration - linear unit
sigPower_noNoise_LR_lin = 1.*sigPower_noNoise_LR_cohe + 1.*sigPower_noNoise_LR_inco;
sigPower_noNoise_RR_lin = 1.*sigPower_noNoise_RR_cohe + 1.*sigPower_noNoise_RR_inco;

%% Test plot - amir - DDM in physical units

% delay_step = geoSYSp.mpas(indexTxFrequency);
% delay_bins = 0 : (geoSYSp.nx_Rdelay - 1);
% DelayAxis_usec = (delay_bins .* delay_step) - geoSYSp.delayOffset_microsec(indexTxFrequency);
% 
% doppler_step = geoSYSp.fpas;
% DopplerAxis_Hz = Freq;
% 
% DDM_Data = sigPower_noNoise_LR_lin;
% 
% grid_X_edges = [DelayAxis_usec - delay_step/2, DelayAxis_usec(end) + delay_step/2];
% grid_Y_edges = [DopplerAxis_Hz - doppler_step/2, DopplerAxis_Hz(end) + doppler_step/2];
% 
% figure('Color','w');
% imagesc(DelayAxis_usec, DopplerAxis_Hz, DDM_Data);
% axis xy;
% colormap parula; colorbar;
% hold on;
% 
% for k = 1:length(grid_X_edges)
%     x_val = grid_X_edges(k);
%     line([x_val x_val], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 0.5);
% end
% 
% for k = 1:length(grid_Y_edges)
%     y_val = grid_Y_edges(k);
%     line(xlim, [y_val y_val], 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 0.5);
% end
% 
% xlabel('Delay [\mus]');
% ylabel('Doppler Frequency [Hz]');
% title('DDM (LR) - Physical Units with Bin Grid');
% set(gca, 'Layer', 'top');
% hold off;

%%

%% Coherent Channel implementation

% indx_origin = surfParam.Xplot == 0.0;
% indy_origin = surfParam.Yplot == 0.0;

%Non ricordo più perchè qui si cerca, il sistema di rif è centrato su %SP
indx_origin = TRAp.indx_origin;
indy_origin = TRAp.indy_origin;

indx_facetAroundSP = surfParam.Xplot < 500 & surfParam.Xplot > -500;
indy_facetAroundSP = surfParam.Xplot < 500 & surfParam.Xplot > -500;

rho_TRsp = DDMparam.TxPmod(indy_facetAroundSP, indx_facetAroundSP) + DDMparam.RxPmod(indy_facetAroundSP, indx_facetAroundSP);
rho_TRsp = mean(rho_TRsp(:));
rho_TR   = DDMparam.TxPmod + DDMparam.RxPmod;

% WAFxy  = S(f-f0)*G(t-t0) 
WAFxy  = 1;
Gt     = TRAp.GainTx_linear(:,:,indexTxFrequency); %%TRAp.TxGainAtSP
Gtmax  = TRAp.TxGainAtSP(indexTxFrequency); %max(Gt(:)); 

factAcoh  = mean(mean(DDMparam.regGridSim_cellareas(indy_facetAroundSP, indx_facetAroundSP))).*TRAp.waveL(indexTxFrequency)*sqrt(Gtmax.*db2pow(TRAp.TxSignalPower(indexTxFrequency)))/sqrt(4*pi); 

factAcoh2 = sqrt(DDMparam.Gain_LHCP_onsurf_lin(indy_facetAroundSP, indx_facetAroundSP,indexTxFrequency).*...
                sigmaLR_cohe(indy_facetAroundSP, indx_facetAroundSP))*WAFxy*...                 
                1./(DDMparam.TxPmod(indy_facetAroundSP, indx_facetAroundSP) + DDMparam.RxPmod(indy_facetAroundSP, indx_facetAroundSP)).*...
                exp(1j*2*pi/TRAp.waveL(indexTxFrequency)*(DDMparam.TxPmod(indy_facetAroundSP, indx_facetAroundSP)  + DDMparam.RxPmod(indy_facetAroundSP, indx_facetAroundSP)));

% AcohLR    = factAcoh.*sum(factAcoh2(:));
AcohLR = factAcoh .* sum(factAcoh2(:), 'omitnan'); %Ignor Nan values Due to low simulation resolution (Grass) - AMIRREZA

factAincoh = DDMparam.regGridSim_cellareas.*TRAp.waveL(indexTxFrequency)^2*Gtmax.*db2pow(TRAp.TxSignalPower(indexTxFrequency))/(4*pi)^3; 

WAFxy_full = 1; %This should be the one on the surface - I do not have it above but sigmaLR_inco is very small.

%evaluate the number of facet within 1 chip - not sure, different unit?
%DDMparam.fdelay is too much variable
nfacet  = floor(geoSYSp.chipTime(indexTxFrequency)/DDMparam.fdelay(indy_origin,indx_origin));
nfaceth = abs(floor(nfacet/2));

nfaceth = 6;

indx_origin2 = find(indx_origin == 1);
indy_origin2 = find(indy_origin == 1);

indx_origin3 = (indx_origin2 - nfaceth):(indx_origin2 + nfaceth);
indy_origin3 = (indy_origin2 - nfaceth):(indy_origin2 + nfaceth);

%it gets only the NRCS around the SP and up to one chip
sigmaLR_inco_1chip = sigmaLR_inco(indy_origin3, indx_origin3);

%AdiffLR0   = factAincoh.*DDMparam.Gain_LHCP_onsurf_lin(:,:,indexTxFrequency).*sigmaLR_inco.*WAFxy_full^2.*...
%                         1./(DDMparam.TxPmod.^2.*DDMparam.RxPmod.^2);
AdiffLR0   = factAincoh(indy_origin3, indx_origin3).*DDMparam.Gain_LHCP_onsurf_lin(indy_origin3, indx_origin3, indexTxFrequency).*sigmaLR_inco_1chip.*WAFxy_full^2.*...
                          1./(DDMparam.TxPmod(indy_origin3, indx_origin3).^2.*DDMparam.RxPmod(indy_origin3, indx_origin3).^2);

%sum over about 1 chip
AdiffLR    = sum(sum(AdiffLR0,'omitnan'));


%% ----------THERMAL AND SPEKLE NOISE------------------------

%generation of N random numbers in the interval (a,b) -> r = a + (b-a).*rand(N,1).

geoSYSp.LN_lin_LHCP = 10.^(instrumentConfParam.LN_dB_LHCP(indexTxFrequency)/10);
geoSYSp.LN_lin_RHCP = 10.^(instrumentConfParam.LN_dB_RHCP(indexTxFrequency)/10);
geoSYSp.LN_lin_RHCP_zenith = 10.^(instrumentConfParam.LN_dB_RHCP_zenith(indexTxFrequency)/10);

geoSYSp.LS_lin_LHCP = 10.^(instrumentConfParam.LS_dB_LHCP(indexTxFrequency)/10);
geoSYSp.LS_lin_RHCP = 10.^(instrumentConfParam.LS_dB_RHCP(indexTxFrequency)/10);
geoSYSp.LS_lin_RHCP_zenith = 10.^(instrumentConfParam.LS_dB_RHCP_zenith(indexTxFrequency)/10);

%physical temperature of the receiver in Kelvin
%geoSYSp.TRx_K         = TRx_K;
DeltaTRx_interf_K = geoSYSp.DeltaTRx_interf_K;

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1 && indexTxFrequency==1
   %if monte carlo is activated it is a random error
   geoSYSp.DeltaTRx_K = -abs(geoSYSp.DeltaTRx_interf_K) + (abs(geoSYSp.DeltaTRx_interf_K) + abs(geoSYSp.DeltaTRx_interf_K))*rand(1);
    elseif tag.MC == 0 && indexTxFrequency==1
   geoSYSp.DeltaTRx_K = geoSYSp.DeltaTRx_interf_K;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
    geoSYSp.DeltaTRx_K = geoSYSp.DeltaTRx_interf_K;
end

%%noise figure LNA - Realization in dB
geoSYSp.NoiseFigure_lin_LHCP        = 10.^(instrumentConfParam.gNF_LHCP(indexTxFrequency)/10);
geoSYSp.NoiseFigure_lin_RHCP        = 10.^(instrumentConfParam.gNF_RHCP(indexTxFrequency)/10);
geoSYSp.NoiseFigure_lin_RHCP_zenith = 10.^(instrumentConfParam.gNF_RHCP_zenith(indexTxFrequency)/10);
geoSYSp.NoiseFigure_dB_LHCP         = instrumentConfParam.gNF_LHCP(indexTxFrequency);
geoSYSp.NoiseFigure_dB_RHCP         = instrumentConfParam.gNF_RHCP(indexTxFrequency);
geoSYSp.NoiseFigure_dB_RHCP_zenith  = instrumentConfParam.gNF_RHCP_zenith(indexTxFrequency);

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
     %if monte carlo is activated it is a random error
     DeltaF_dB        = -abs(geoSYSp.DeltaF_interf_dB) + (abs(geoSYSp.DeltaF_interf_dB) + abs(geoSYSp.DeltaF_interf_dB))*rand(1);
    else
        DeltaF_dB        = geoSYSp.DeltaF_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaF_dB        = geoSYSp.DeltaF_interf_dB;
end

geoSYSp.DeltaF_lin   = 10.^(DeltaF_dB/10);
geoSYSp.DeltaF_dB    = DeltaF_dB;

%system gain
geoSYSp.G_sys0_lin_LHCP      = 10.^(instrumentConfParam.Gsys_dB_LHCP(indexTxFrequency)/10);
geoSYSp.G_sys0_lin_RHCP      = 10.^(instrumentConfParam.Gsys_dB_RHCP(indexTxFrequency)/10);
geoSYSp.G_sys0_lin_RHCP_zenith      = 10.^(instrumentConfParam.Gsys_dB_RHCP_zenith(indexTxFrequency)/10);
geoSYSp.G_sys0_dB_LHCP       = instrumentConfParam.Gsys_dB_LHCP(indexTxFrequency);
geoSYSp.G_sys0_dB_RHCP       = instrumentConfParam.Gsys_dB_RHCP(indexTxFrequency);
geoSYSp.G_sys0_dB_RHCP_zenith       = instrumentConfParam.Gsys_dB_RHCP_zenith(indexTxFrequency);

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaG_sys_dB      = -abs(geoSYSp.DeltaG_sys_interf_dB) + (abs(geoSYSp.DeltaG_sys_interf_dB) + abs(geoSYSp.DeltaG_sys_interf_dB))*rand(1);
    else
        DeltaG_sys_dB      = geoSYSp.DeltaG_sys_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaG_sys_dB      = geoSYSp.DeltaG_sys_interf_dB;
end

geoSYSp.DeltaG_sys_lin = 10.^(DeltaG_sys_dB/10);
geoSYSp.DeltaG_sys_dB  = DeltaG_sys_dB;

%antenna temperature and receiver temperature
geoSYSp.ReceiverTemperature_K = instrumentConfParam.Tref_K;

%coefficients in linear units
% coefficients Noise Figure
%geoSYSp.gNF = instrumentConfParam.gNF(indexTxFrequency);
geoSYSp.mNF_LHCP = instrumentConfParam.mNF_LHCP(indexTxFrequency);
geoSYSp.mNF_RHCP = instrumentConfParam.mNF_RHCP(indexTxFrequency);
geoSYSp.mNF_RHCP_zenith = instrumentConfParam.mNF_RHCP_zenith(indexTxFrequency);

% coefficients Gain
geoSYSp.gG_LHCP = instrumentConfParam.gG_LHCP(indexTxFrequency);
geoSYSp.gG_RHCP = instrumentConfParam.gG_RHCP(indexTxFrequency);
%geoSYSp.gG_RHCP_zenith = instrumentConfParam.gG_RHCP_zenith(indexTxFrequency);

geoSYSp.mG_LHCP = instrumentConfParam.mG_LHCP(indexTxFrequency);
geoSYSp.mG_RHCP = instrumentConfParam.mG_RHCP(indexTxFrequency);
geoSYSp.mG_RHCP_zenith = instrumentConfParam.mG_RHCP_zenith(indexTxFrequency);


%% fluctuations noise bb, physical temperature, and gain of the receiver

%%%Receiver temperature
%realization based on the same interval

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1 && indexTxFrequency==1
        %if monte carlo is activated it is a random error
        geoSYSp.DeltaTRxbb_K         = -abs(DeltaTRx_interf_K)+ (abs(DeltaTRx_interf_K)+ abs(DeltaTRx_interf_K))*rand(1);
    elseif tag.MC == 0 && indexTxFrequency==1
        geoSYSp.DeltaTRxbb_K         = DeltaTRx_interf_K;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
geoSYSp.DeltaTRxbb_K         = DeltaTRx_interf_K;
end

%geoSYSp.DeltaTRxbb_K  = DeltaTRxbb_K;

%%%Load temperature in K, same error of the Receiver temperature
%realization based on the same interval

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1 && indexTxFrequency==1
        %if monte carlo is activated it is a random error
        geoSYSp.DeltaTL_K         = -abs(DeltaTRx_interf_K) + (abs(DeltaTRx_interf_K) + abs(DeltaTRx_interf_K))*rand(1);
    elseif tag.MC == 0 && indexTxFrequency==1
        geoSYSp.DeltaTL_K         = DeltaTRx_interf_K;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
geoSYSp.DeltaTL_K         = DeltaTRx_interf_K;
end

%geoSYSp.DeltaTL_K     = DeltaTL_K;

%realization based on the same interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaFbb             = -abs(geoSYSp.DeltaF_interf_dB)     + (abs(geoSYSp.DeltaF_interf_dB)     + abs(geoSYSp.DeltaF_interf_dB))*rand(1);
    else
        DeltaFbb             = geoSYSp.DeltaF_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaFbb             = geoSYSp.DeltaF_interf_dB;
end

geoSYSp.DeltaFbb_lin     = 10.^(DeltaFbb/10);
geoSYSp.DeltaFbb_dB      = DeltaFbb;

%realization based on the same interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaG_sysbb         = -abs(geoSYSp.DeltaG_sys_interf_dB) + (abs(geoSYSp.DeltaG_sys_interf_dB) + abs(geoSYSp.DeltaG_sys_interf_dB))*rand(1);
    else
        DeltaG_sysbb         = geoSYSp.DeltaG_sys_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaG_sysbb         = geoSYSp.DeltaG_sys_interf_dB;
end

geoSYSp.DeltaG_sysbb_lin = 10.^(DeltaG_sysbb/10);
geoSYSp.DeltaG_sysbb_dB  = DeltaG_sysbb;

%%%%System and noise characterization
[sigPower_Noise_LR_lin, sigPower_Noise_RR_lin, NLbb_LR, NLbb_RR, G_sys_lin_LHCP, G_sys_lin_RHCP, G_sys_lin_LHCP_nominal , G_sys_lin_RHCP_nominal] = gen_noise_sys(sigPower_noNoise_LR_lin, sigPower_noNoise_RR_lin, geoSYSp);
geoSYSp.G_sys_lin_LHCP(indexTxFrequency) = G_sys_lin_LHCP;
geoSYSp.G_sys_lin_RHCP(indexTxFrequency) = G_sys_lin_RHCP;
geoSYSp.G_sys_lin_LHCP_nominal(indexTxFrequency) = G_sys_lin_LHCP_nominal;
geoSYSp.G_sys_lin_RHCP_nominal(indexTxFrequency) = G_sys_lin_RHCP_nominal;

%% ----------- UP LOOKING ANTENNA POWER -----------------------------------
%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaTRx_up                = -abs(DeltaTRx_interf_K)   + (abs(DeltaTRx_interf_K) + abs(DeltaTRx_interf_K))*rand(1);
    else
        DeltaTRx_up                = DeltaTRx_interf_K;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaTRx_up                = DeltaTRx_interf_K;
end

geoSYSp.DeltaTRx_up_K          = DeltaTRx_up;

%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaG_sys_up             = -abs(geoSYSp.DeltaG_sys_interf_dB) + (abs(geoSYSp.DeltaG_sys_interf_dB) + abs(geoSYSp.DeltaG_sys_interf_dB))*rand(1);
    else
        DeltaG_sys_up             = geoSYSp.DeltaG_sys_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaG_sys_up             = geoSYSp.DeltaG_sys_interf_dB;
end

geoSYSp.DeltaG_sys_up_lin      = 10.^(DeltaG_sys_up/10);
geoSYSp.DeltaG_sys_up_dB       = DeltaG_sys_up;

%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaF_up                  = -abs(geoSYSp.DeltaF_interf_dB) + (abs(geoSYSp.DeltaF_interf_dB) + abs(geoSYSp.DeltaF_interf_dB))*rand(1);
    else
        DeltaF_up                  = geoSYSp.DeltaF_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaF_up                  = geoSYSp.DeltaF_interf_dB;
end

geoSYSp.DeltaF_up_lin          = 10.^(DeltaF_up/10);
geoSYSp.DeltaF_up_dB           = DeltaF_up;

%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaTRxbb_up              = -abs(DeltaTRx_interf_K) + (abs(DeltaTRx_interf_K) + abs(DeltaTRx_interf_K))*rand(1);
    else
        DeltaTRxbb_up              = DeltaTRx_interf_K;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaTRxbb_up              = DeltaTRx_interf_K;
end
geoSYSp.DeltaTRxbb_up          = DeltaTRxbb_up;
 
%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaFbb_up                = -abs(geoSYSp.DeltaF_interf_dB)      + (abs(geoSYSp.DeltaF_interf_dB) + abs(geoSYSp.DeltaF_interf_dB))*rand(1);
    else
        DeltaFbb_up                = geoSYSp.DeltaF_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaFbb_up                = geoSYSp.DeltaF_interf_dB;
end

geoSYSp.DeltaFbb_up_lin        = 10.^(DeltaFbb_up/10);
geoSYSp.DeltaFbb_up_dB         = DeltaFbb_up;

%realization within the prescribed interval
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    if tag.MC == 1
        %if monte carlo is activated it is a random error
        DeltaG_sysbb_up            = -abs(geoSYSp.DeltaG_sys_interf_dB) + (abs(geoSYSp.DeltaG_sys_interf_dB) + abs(geoSYSp.DeltaG_sys_interf_dB))*rand(1);
    else
        DeltaG_sysbb_up            = geoSYSp.DeltaG_sys_interf_dB;
    end
elseif isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
DeltaG_sysbb_up            = geoSYSp.DeltaG_sys_interf_dB;
end

geoSYSp.DeltaG_sysbb_up_lin    = 10.^(DeltaG_sysbb_up/10);
geoSYSp.DeltaG_sysbb_up_dB     = DeltaG_sysbb_up;
 

[LR_UP_Power, RR_UP_Power, LR_UP_Power_Flu, RR_UP_Power_Flu, PN_star, NLbb_UP, UP_pattern, TxGainAtRx] = UP_power(TRAp, geoSYSp, indexTxFrequency);


%% ----------------- DDM plots -----------------------------------------

Fplot = 1:1:geoSYSp.ny_Dfreq;
Rplot = 1:1:geoSYSp.nx_Rdelay;
% Rplot = Rplot * geoSYSp.mpas(1); %Amir
% Fplot = Fplot * geoSYSp.fpas; %amir
%Is this really needed?
sigPower_noNoise_LR_lin(sigPower_noNoise_LR_lin <= 0) = 0; 
sigPower_Noise_LR_lin(sigPower_Noise_LR_lin     <= 0) = 0;

%db converting
sigPower_noNoise_LR_dB = db(sigPower_noNoise_LR_lin,'power');
sigPower_Noise_LR_dB   = db(sigPower_Noise_LR_lin,'power');

%max computing
max_sigPower_noNoise_LR = max(max(sigPower_noNoise_LR_lin));
max_sigPower_Noise_LR = max(max(sigPower_Noise_LR_lin));

%test power
%10*log10(max_sigPower_noNoise_LR)

% %with test antenna small numbers can introduce a non-zero imaginary part
% if sum(imag(sigPower_noNoise_LR_dB(:))) ~= 0
%     sigPower_noNoise_LR_dB  = real(sigPower_noNoise_LR_dB);
%     disp('output LR signal power has complex values, please check'),
% end

% plot of power DDM, noise-free and noisy
if tag.plotDDM == 1

    if isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
    figDDM_LRtotPower = figure('Name', 'LR power clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,sigPower_noNoise_LR_dB);
    colorbar
    colormap('jet')
    axis xy
    caxis([max(sigPower_noNoise_LR_dB(:)) - 20, max(sigPower_noNoise_LR_dB(:))])

    if indexTxFrequency==1
    title(' LR1 power clean DDM in dB Watt ')
    elseif indexTxFrequency==2
    title(' LR5 power clean DDM in dB Watt ')
    end

    figDDM_LRtotPower_lin = figure('Name', 'LR power clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(sigPower_noNoise_LR_dB));
    colorbar
    colormap('parula')
    axis xy

    if indexTxFrequency==1
    title(' LR1 power clean DDM in Watt ')
    elseif indexTxFrequency==2
    title(' LR5 power clean DDM in Watt ')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')
    caxis([0, max(db2pow(sigPower_noNoise_LR_dB(:)))])

    if indexTxFrequency==1
        saveas(figDDM_LRtotPower, strjoin([nameOutputFile '_DDMnoiseFree_LR1.fig'],''));
        saveas(figDDM_LRtotPower_lin, strjoin([nameOutputFile '_DDMnoiseFree_LR1_lin.fig'],''));
    else
        saveas(figDDM_LRtotPower, strjoin([nameOutputFile '_DDMnoiseFree_LR5.fig'],''));
        saveas(figDDM_LRtotPower_lin, strjoin([nameOutputFile '_DDMnoiseFree_LR5_lin.fig'],''));
    end

    end


    if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    figDDM_LRtotPowerNoisy = figure('Name', 'LR power noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,sigPower_Noise_LR_dB);
    colorbar
    colormap('jet')
    axis xy
    caxis([max(sigPower_Noise_LR_dB(:)) - 20, max(sigPower_Noise_LR_dB(:))])

    if indexTxFrequency==1
    title(' LR1 power noisy DDM in dB Watt ')
    elseif indexTxFrequency==2
    title(' LR5 power noisy DDM in dB Watt ')
    end
    
    figDDM_LRtotPowerNoisy_lin = figure('Name', 'LR power noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(sigPower_Noise_LR_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(sigPower_Noise_LR_dB(:)))])

    if indexTxFrequency==1
    title(' LR1 power noisy DDM in Watt ')
    elseif indexTxFrequency==2
    title(' LR5 power noisy DDM in Watt ')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')

    if indexTxFrequency==1
        saveas(figDDM_LRtotPowerNoisy, strjoin([nameOutputFile '_DDMnoisy_LR1.fig'],''));
        saveas(figDDM_LRtotPowerNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_LR1_lin.fig'],''));
    else
        saveas(figDDM_LRtotPowerNoisy, strjoin([nameOutputFile '_DDMnoisy_LR5.fig'],''));
        saveas(figDDM_LRtotPowerNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_LR5_lin.fig'],''));
    end
    
    end


end


%%%is this needed?
sigPower_noNoise_RR_lin(sigPower_noNoise_RR_lin <= 0) = 0; 
sigPower_Noise_RR_lin(sigPower_Noise_RR_lin <= 0)     = 0;

%dB converting
sigPower_noNoise_RR_dB = db(sigPower_noNoise_RR_lin,'power');
sigPower_Noise_RR_dB   = db(sigPower_Noise_RR_lin,'power');

%max computing
max_sigPower_noNoise_RR  = max(max(sigPower_noNoise_RR_lin));
max_sigPower_Noise_RR    = max(max(sigPower_Noise_RR_lin));

% %with test antenna small numbers can introduce a non-zero imaginary part
% if sum(imag(sigPower_noNoise_RR_dB(:))) ~= 0
%     sigPower_noNoise_RR_dB  = real(sigPower_noNoise_RR_dB);
%     disp('output LR signal power has complex values, please check'),
% end

% plot of power DDM, noise-free and noisy
if tag.plotDDM == 1

    if isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
    figDDM_RRtotPower = figure('Name', 'RR power clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,sigPower_noNoise_RR_dB);
    colorbar
    colormap('jet')
    axis xy
    caxis([max(sigPower_noNoise_RR_dB(:)) - 20, max(sigPower_noNoise_RR_dB(:))])

    if indexTxFrequency==1
    title(' RR1 power clean DDM in dB Watt ')
    elseif indexTxFrequency==2
    title(' RR5 power clean DDM in dB Watt ')
    end

    figDDM_RRtotPower_lin = figure('Name', 'RR power clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(sigPower_noNoise_RR_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(sigPower_noNoise_RR_dB(:)))])

    if indexTxFrequency==1
    title(' RR1 power clean DDM in Watt ')
    elseif indexTxFrequency==2
    title(' RR5 power clean DDM in Watt ')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')

    if indexTxFrequency==1
        saveas(figDDM_RRtotPower, strjoin([nameOutputFile '_DDMnoiseFree_RR1.fig'],''));
        saveas(figDDM_RRtotPower_lin, strjoin([nameOutputFile '_DDMnoiseFree_RR1_lin.fig'],''));
    else
        saveas(figDDM_RRtotPower, strjoin([nameOutputFile '_DDMnoiseFree_RR5.fig'],''));
        saveas(figDDM_RRtotPower_lin, strjoin([nameOutputFile '_DDMnoiseFree_RR5_lin.fig'],''));
    end

    end
    
    if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    figDDM_RRtotPowerNoisy = figure('Name', 'RR power noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,sigPower_Noise_RR_dB);
    colorbar
    colormap('jet')
    axis xy
    caxis([max(sigPower_Noise_RR_dB(:)) - 20, max(sigPower_Noise_RR_dB(:))])

    if indexTxFrequency==1
    title(' RR1 power noisy DDM in dB Watt ')
    elseif indexTxFrequency==2
    title(' RR5 power noisy DDM in dB Watt ')
    end

    figDDM_RRtotPowerNoisy_lin = figure('Name', 'RR power noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(sigPower_Noise_RR_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(sigPower_Noise_RR_dB(:)))])

    xlabel('Delay bin')
    ylabel('Doppler bin')
    
    if indexTxFrequency==1
    title(' RR1 power noisy DDM in Watt ')
    elseif indexTxFrequency==2
    title(' RR5 power noisy DDM in Watt ')
    end

    if indexTxFrequency==1
        saveas(figDDM_RRtotPowerNoisy, strjoin([nameOutputFile '_DDMnoisy_RR1.fig'],''));
        saveas(figDDM_RRtotPowerNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_RR1_lin.fig'],''));
    else
        saveas(figDDM_RRtotPowerNoisy, strjoin([nameOutputFile '_DDMnoisy_RR5.fig'],''));
        saveas(figDDM_RRtotPowerNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_RR5_lin.fig'],''));
    end

    end

end


%% %%%% COVARIANCE MODULE

if tag.Cov == 1
    addpath('DDM_Covariance_Speckle_and_ThermalNoise_Matlab_Version')
    Main_Generate_DDM_Cov
end

end

