%% ************************************************************************
% MODULE NAME:      ionoDelayComputation.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, D. Comite, N. Pierdicca, L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      Oct 2022 by Laura Dente
% Main updates:          Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module prepares all the inputs to call the routine runNTCM, to
% compute the ionospheric delay alonf Tx-SP and Rx-SP path
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function ionosphericDelay = ionoDelayComputation(ionisationLevelParam, GNSS_Week, GNSS_SoW,...
                                                RxECEF,TxECEF,SP_lla,TRAp)

%% ----- prepare the input matrix of the NTCM model -----------------------

% convert GNSS_Week and GNSS_Sow in DOY and UTC hour
totalGNSStime = datenum(1980, 1, 6, 0 , 0 ,0)+GNSS_Week*7+GNSS_SoW/86400 ;
UTCtime=datetime(totalGNSStime, 'ConvertFrom', 'datenum') ;
GNSS_doy=day(UTCtime,'dayofyear');
GNSS_hour=hour(UTCtime);

% convert ECEF coordinates to geodetic coordinates
RX_lla = ecef2lla(RxECEF); %(1,3) latitude, longitude and altitude
TX_lla = ecef2lla(TxECEF); %(1,3) latitude, longitude and altitude

% signal frequency (Hz)
carrFreq=TRAp.TxFrequency.*(10^9);

% set input matrix to compute the ionospheric delay along Tx-SP path
inputTxSP(1,1) = ionisationLevelParam(1); % ai0 medium solar activity
inputTxSP(1,2) = ionisationLevelParam(2); % ai1 medium solar activity
inputTxSP(1,3) = ionisationLevelParam(3); % ai2 medium solar activity
inputTxSP(1,4) = GNSS_doy; % DOY
inputTxSP(1,5) = GNSS_hour; % UTC hour
inputTxSP(1,6) = SP_lla(1,2); % User receiver Geodetic longitude
inputTxSP(1,7) = SP_lla(1,1); % User receiver Geodetic latitude
inputTxSP(1,8) = SP_lla(1,3); % User receiver Geodetic height
inputTxSP(1,9) = TX_lla(1,2); % Satellite Geodetic longitude (TX)
inputTxSP(1,10)= TX_lla(1,1); % Satellite Geodetic latitude (TX)
inputTxSP(1,11)= TX_lla(1,3); % Satellite Geodetic height Range (TX)

% set input matrix to compute the ionospheric delay along Tx-SP path
inputRxSP(1,1) = ionisationLevelParam(1); % ai0 medium solar activity
inputRxSP(1,2) = ionisationLevelParam(2); % ai1 medium solar activity
inputRxSP(1,3) = ionisationLevelParam(3); % ai2 medium solar activity
inputRxSP(1,4) = GNSS_doy; % DOY
inputRxSP(1,5) = GNSS_hour; % UTC hour
inputRxSP(1,6) = SP_lla(1,2); % User receiver Geodetic longitude
inputRxSP(1,7) = SP_lla(1,1); % User receiver Geodetic latitude
inputRxSP(1,8) = SP_lla(1,3); % User receiver Geodetic height
inputRxSP(1,9) = RX_lla(1,2); % Satellite Geodetic longitude (RX)
inputRxSP(1,10)= RX_lla(1,1); % Satellite Geodetic latitude (RX)
inputRxSP(1,11)= RX_lla(1,3); % Satellite Geodetic height Range (RX)

for indexTxFrequency=1:length(carrFreq)
    %% ------ run the NTCM model for the path Tx-SP ------
    [~, ~, IonoDelayTxSP_m] = runNTCM(inputTxSP, carrFreq(indexTxFrequency)); %meters

    %% ------ run the NTCM model for the path Rx-SP ------
    [~, ~, IonoDelayRxSP_m] = runNTCM(inputRxSP, carrFreq(indexTxFrequency)); %meters

    %% ------ compute the total delay for the path Tx-SP-Rx
%     ionosphericDelay(indexTxFrequency) = IonoDelayTxSP_m;+IonoDelayRxSP_m; % Original

% [AMIRREZA]
% In the Great project, due to the low altitude of the airborne platform,
% there is no ionospheric delay for the receiver (Rx).
% In the original code, the function returned NaN for Rx.
% use omitnan to avoid this issue

    ionosphericDelay(indexTxFrequency) = sum([IonoDelayTxSP_m, IonoDelayRxSP_m], 'omitnan');

end

end