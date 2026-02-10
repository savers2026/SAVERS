%% ************************************************************************
% MODULE NAME:      HSAVERS_IHS_OENUonDEM.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, D. Comite, N. Pierdicca, L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      2010 by Roberto Giusto (airborne)
% Main updates:          2018 – 2022 by L. Dente (spaceborne, inhomogenoues cover, HSAVERS)
%                        2021 - 2022 by D. Comite (optimization, HSAVERS)
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% This module calls the routines to compute the specular point position, 
% to read the DEM and to land cover map, to compute the ionospheric delay,
% to compute the parameters at each facet,
% to run the Fortran routien (AIEM and forest model), to compute the Tx gain at SP
% and to compute the DDMs.
% It transforms the reference system from ECEF to SAVERS Global Refence system
% (rotated ENU), it computes the reflectivit and it evaluates the coherent
% component.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************n

function [outputSimulations, infoAtSPforRefFile, UP_pattern] = HSAVERS_IHS_OENUonDEM(dest_dir, namefile, TxECEF, VTxECEF, RxECEF, VRxECEF,...
                                                                   geoSYSp, TRAp, DEMinfo, landCoverMap, coverMap_IDs, instrumentConfParam, terrainParameters,...
                                                                   vegetationParameters, ionisationLevelParam, GNSS_Week, GNSS_SoW, tag,landCoverMap_all)



%goes in debug mod only if not generating exe
if tag.exe == 0, dbstop if error, end

%Delete existing files
fclose('all');

%% ------- READ Tx and Rx POSITION AND VELOCITY IN ECEF -------------------
%count the number of Rx positions along the selected track (= no of specular points)
numberOfRxPos = size(RxECEF);
%The simulation will be repeated Nruns times, i.e. for all SP
Nruns         = numberOfRxPos(1);

%% ------->  Specular Point Calculation over WGS84 (no elevation)  <-------
 
%find the SP position using the method of the minimum distance Tx-SP-Rx
%the computation is carried out assuming a flat surface on WGS84 ellipsoid

%Initializing variables
SP_lla_onEllipsoid = zeros(Nruns,3);

for i = 1:Nruns         % for each SP
    [SP_lla_onEllipsoid_i, ~] = Specular_Point_Position_Calculation(RxECEF(i,:), TxECEF(i,:));

    %----\LLH coordinates of Specular Point on the WGS84 ellipsoid
    SP_lla_onEllipsoid(i,:)   = SP_lla_onEllipsoid_i(:);
    
%     %Evaluate difference between info in Orbit file and the one computed here
%     SPlat_diff(i) = geoSYSp.SPlat_orbitFile(i) - SP_lla_onEllipsoid_i(1);
%     Splon_diff(i) = geoSYSp.SPlon_orbitFile(1) - SP_lla_onEllipsoid_i(2);    
end

%% Pre-allocation of variables
%Defining specular point variables in LLA and ECEF
SP_lla           = zeros(Nruns,3);
%SP_ECEF_ElevZero = zeros(Nruns,3);
SP_ECEF_meanElev = zeros(Nruns,3);
SPelevFromDEM = zeros(Nruns,3);

max_sigPower_noNoise_LR = zeros(Nruns,2);
max_sigPower_noNoise_RR = zeros(Nruns,2);
max_sigPower_Noise_LR   = zeros(Nruns,2);
max_sigPower_Noise_RR   = zeros(Nruns,2);
max_reflectivity_noNoise_LR_tot_dB = zeros(Nruns,2);
max_reflectivity_noNoise_RR_tot_dB = zeros(Nruns,2);
max_reflectivity_Noisy_LR_tot_dB = zeros(Nruns,2);
max_reflectivity_Noisy_RR_tot_dB = zeros(Nruns,2);
max_geobReflec_noNoise_LR_dB = zeros(Nruns,2);
max_geobReflec_noNoise_RR_dB = zeros(Nruns,2);
max_geobPower_Noisy_LR = zeros(Nruns,2);
max_geobPower_Noisy_RR = zeros(Nruns,2);

sigPower_noNoise_LR_lin = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2); 
sigPower_noNoise_RR_lin = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2); 
sigPower_Noise_LR_lin   = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
sigPower_Noise_RR_lin   = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);

incAngle_enuRot_vect = zeros(1,Nruns);

%UP-looking Antenna Power          
LR_UP_Power = zeros(Nruns,2);
RR_UP_Power = zeros(Nruns,2);
UP_pattern  = zeros(Nruns,2);

LR_UP_Power_Flu = zeros(Nruns,2);
RR_UP_Power_Flu = zeros(Nruns,2);
PN_star         = zeros(Nruns,2);

NLbb_DOWN_LR    = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
NLbb_DOWN_RR    = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
NLbb_UP         = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);

SNR_LR       = zeros(Nruns,2);
SNR_RR       = zeros(Nruns,2);
geobSNR_LR   = zeros(Nruns,2);
geobSNR_RR   = zeros(Nruns,2);

delayAtSPReferencedToDirect_1_microsec = zeros(Nruns,1);
delayAtSPReferencedToDirect_5_microsec = zeros(Nruns,1);
absoluteDelayAtSP_1_microsec           = zeros(Nruns,1);
absoluteDelayAtSP_5_microsec           = zeros(Nruns,1);

%% modified for MM & Ideal antenna simultaneously (Amir)

max_sigPower_noNoise_NoError_LR = zeros(Nruns,2);
max_sigPower_noNoise_NoError_RR = zeros(Nruns,2);
max_sigPower_Noise_NoError_LR   = zeros(Nruns,2);
max_sigPower_Noise_NoError_RR   = zeros(Nruns,2);
max_reflectivity_noNoise_NoError_LR_tot_dB = zeros(Nruns,2);
max_reflectivity_noNoise_NoError_RR_tot_dB = zeros(Nruns,2);
max_reflectivity_Noisy_NoError_LR_tot_dB = zeros(Nruns,2);
max_reflectivity_Noisy_NoError_RR_tot_dB = zeros(Nruns,2);
max_geobReflec_noNoise_NoError_LR_dB = zeros(Nruns,2);
max_geobReflec_noNoise_NoError_RR_dB = zeros(Nruns,2);
max_geobPower_Noisy_NoError_LR = zeros(Nruns,2);
max_geobPower_Noisy_NoError_RR = zeros(Nruns,2);
sigPower_noNoise_NoError_LR_lin = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2); 
sigPower_noNoise_NoError_RR_lin = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2); 
sigPower_Noise_NoError_LR_lin   = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
sigPower_Noise_NoError_RR_lin   = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);

LR_UP_Power_NoError = zeros(Nruns,2);
RR_UP_Power_NoError = zeros(Nruns,2);
UP_pattern_NoError  = zeros(Nruns,2);

LR_UP_Power_Flu_NoError = zeros(Nruns,2);
RR_UP_Power_Flu_NoError = zeros(Nruns,2);

NLbb_DOWN_NoError_LR    = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
NLbb_DOWN_NoError_RR    = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);
NLbb_UP_NoError         = zeros(geoSYSp.ny_Dfreq, geoSYSp.nx_Rdelay, Nruns, 2);

SNR_NoError_LR       = zeros(Nruns,2);
SNR_NoError_RR       = zeros(Nruns,2);
PN_star_NoError         = zeros(Nruns,2);


%define a flooding flag
% if in the land cover class of the site reports the classes 5 an 6
% (flooded forest and flooded low vegetation) and the user has set a flooded
% condition for the surface then the flag is 1
indexFloodedClasses = find(coverMap_IDs>=5);
if indexFloodedClasses & ~isempty(find(terrainParameters(indexFloodedClasses,1)==2))
    flagFlooding = 1;
else
    flagFlooding = 0;
end


%% ---- LOOP over the Nruns, number of SP

for i = 1:Nruns         

    %initializing the variable: nedeed since different runs may have
    %different sizes
    dataDEM_enuRot = [];

    disp(['SP ', num2str(i) , ' of ' , num2str(Nruns)])

    %----------------------------------------------------------
    nameOutputFile = [dest_dir char(namefile(i))];
    %----------------------------------------------------------


    %% No error variables in TRAp (Amir)
    TRAp_NoError = TRAp;
    TRAp_NoError.TxAttitude_angles = [0,0,0];
    TRAp_NoError.RxAttitude_angles = [0,0,0];
    TRAp_NoError.InternalFlag_NoError = 1;

    TRAp_NoError.flagTxMismatch = 0;
    TRAp_NoError.flagRxDOWNmismatch = 0;
    TRAp_NoError.flagRxUPmismatch = 0;

    %% No error variables in geySYSp (Amir)
    geoSYSp_NoError = geoSYSp;
    geoSYSp_NoError.DeltaF_interf_dB = 0;
    geoSYSp_NoError.DeltaTRx_interf_K = 0;
    geoSYSp_NoError.DeltaG_sys_interf_dB = 0;
%     geoSYSp_NoError.alfa = 1;

    %% No error variables in instrumentConfParam (Amir)
    instrumentConfParam_NoError = instrumentConfParam;
%     instrumentConfParam_NoError.gNF_LHCP = [0,0];
%     instrumentConfParam_NoError.gNF_RHCP = [0,0];
%     instrumentConfParam_NoError.gNF_RHCP_zenith = [0,0];
%     instrumentConfParam_NoError.LN_dB_LHCP = [0,0];
%     instrumentConfParam_NoError.LN_dB_RHCP = [0,0];
%     instrumentConfParam_NoError.LN_dB_RHCP_zenith = [0,0];
%     instrumentConfParam_NoError.LS_dB_LHCP = [0,0];
%     instrumentConfParam_NoError.LS_dB_RHCP = [0,0];
%     instrumentConfParam_NoError.LS_dB_RHCP_zenith = [0,0];
%     instrumentConfParam_NoError.mNF_LHCP = [0,0];
%     instrumentConfParam_NoError.mNF_RHCP = [0,0];
%     instrumentConfParam_NoError.mNF_RHCP_zenith = [0,0];
    %% ---------

    %% Definition of the Rx attitude angles vectors

    %TRAp.TxAttitude_angles = [Tx_phi_roll, Tx_theta_pitch, Tx_psi_yaw];

    %if monte carlo is activated the value is made random and overwritten
    if strcmp(geoSYSp.Mode ,'Monte Carlo') == 1
        TRAp.TxAttitude_angles(1) = -abs(TRAp.TxAttitude_angles_GUI(1)) + (abs(TRAp.TxAttitude_angles_GUI(1)) + abs(TRAp.TxAttitude_angles_GUI(1)))*rand(1);
        TRAp.TxAttitude_angles(2) = -abs(TRAp.TxAttitude_angles_GUI(2)) + (abs(TRAp.TxAttitude_angles_GUI(2)) + abs(TRAp.TxAttitude_angles_GUI(2)))*rand(1);
        TRAp.TxAttitude_angles(3) = -abs(TRAp.TxAttitude_angles_GUI(3)) + (abs(TRAp.TxAttitude_angles_GUI(3)) + abs(TRAp.TxAttitude_angles_GUI(3)))*rand(1); 
    else
        TRAp.TxAttitude_angles(1) = TRAp.TxAttitude_angles_GUI(1);
        TRAp.TxAttitude_angles(2) = TRAp.TxAttitude_angles_GUI(2);
        TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles_GUI(3);       
    end

    %if monte carlo is activated the value is made random and overwritten
    if strcmp(geoSYSp.Mode ,'Monte Carlo') == 1
        TRAp.RxAttitude_angles(1) = -abs(TRAp.RxAttitude_angles_GUI(1)) + (abs(TRAp.RxAttitude_angles_GUI(1)) + abs(TRAp.RxAttitude_angles_GUI(1)))*rand(1);
        TRAp.RxAttitude_angles(2) = -abs(TRAp.RxAttitude_angles_GUI(2)) + (abs(TRAp.RxAttitude_angles_GUI(2)) + abs(TRAp.RxAttitude_angles_GUI(2)))*rand(1);
        TRAp.RxAttitude_angles(3) = -abs(TRAp.RxAttitude_angles_GUI(3)) + (abs(TRAp.RxAttitude_angles_GUI(3)) + abs(TRAp.RxAttitude_angles_GUI(3)))*rand(1);
    else
        TRAp.RxAttitude_angles(1) = TRAp.RxAttitude_angles_GUI(1);
        TRAp.RxAttitude_angles(2) = TRAp.RxAttitude_angles_GUI(2);
        TRAp.RxAttitude_angles(3) = TRAp.RxAttitude_angles_GUI(3);
    end


    %% Add Tx power error to Tx power
    for indexTxFreq=1:length(TRAp.TxPower_dBW)
        %% No error variables in TRAp (Amir)
        TRAp_NoError.TxSignalPower(indexTxFreq) = TRAp.TxPower_dBW(indexTxFreq);

        if tag.MC == 1
            %if monte carlo is activated it is a random error %Realization of the fluctuation
            DeltaTxpower_dB = -abs(TRAp.TxPowerError_dBW) + (abs(TRAp.TxPowerError_dBW) + abs(TRAp.TxPowerError_dBW))*rand(1);
        else
            DeltaTxpower_dB = TRAp.TxPowerError_dBW;
        end

        TRAp.TxSignalPower(indexTxFreq)  = TRAp.TxPower_dBW(indexTxFreq) + DeltaTxpower_dB;
    end

  %AMIR  
    %% READ DEM (.MAT OR GEOTIFF FORMAT)        
     [coordDEM_ecef, DEM_lla, resolutionLatDEM, resolutionLonDEM] = readDEM(DEMinfo, SP_lla_onEllipsoid(i,:), flagFlooding);

    % read the 5 km DEM to get the elevation at the SP_lla_onEllipsoid
    if DEMinfo.flatsurface==0
        meanSurfaceElevation = readDEM_5km(DEMinfo, SP_lla_onEllipsoid(i,:));
    else
        meanSurfaceElevation = 0.0;
    end


    %% -------> Specular Point Calculation over an area with average surface elevation  <----------

    %wgs84 = wgs84Ellipsoid('meters');

    %find the SP position using the method of the minimum distance Tx-SP-Rx
    %the computation is carried out assuming a flat surface on a larger WGS84
    %ellipsoid (the axes length is increased by a quantity equal to the mean
    %surface elevation of the area around the SP previuosly found)
    [SP_lla_i,SP_ECEF_i] = Specular_Point_Position_Calculation_PlusMeanElevation(RxECEF(i,:), TxECEF(i,:), meanSurfaceElevation);

    %clear meanSurfaceElevation
    
    SP_lla(i,1) = SP_lla_i(1); % Lat of Specular Point
    SP_lla(i,2) = SP_lla_i(2); % Lon of Specular Point
    SP_lla(i,3) = SP_lla_i(3); % Mean elevation in a window of 10kmx10km around Specular Point 

    SP_ECEF_meanElev(i,1) = SP_ECEF_i(1);
    SP_ECEF_meanElev(i,2) = SP_ECEF_i(2);
    SP_ECEF_meanElev(i,3) = SP_ECEF_i(3);

    clear SP_lla_i
    
    %look for the DEM elevation at the SP coordinates, with the orginal DEM
    %resolution.This info is not needed in HSAVERS, but it is reported in the output
%     indexSPlatInDEM = find((abs(DEM_lla(:,1,1) - SP_lla(i,1))) < resolutionLatDEM/2);
%     indexSPlonInDEM = find((abs(DEM_lla(1,:,2) - SP_lla(i,2))) < resolutionLonDEM/2);
%     SPelevFromDEM(i) = double(DEM_lla(indexSPlatInDEM,indexSPlonInDEM,3)); % Elevation of Specular Point over WGS84 ellipsoid
    indexSPlatInDEM = 1;
    indexSPlonInDEM = 1;
    SPelevFromDEM(i) = 1; % Elevation of Specular Point over WGS84 ellipsoid

    clear indexSPlatInDEM indexSPlonInDEM    

    %projection of the specular point on the ellipsoid, in ECEF coordinates
    %[SP_ECEF_ElevZero(i,1), SP_ECEF_ElevZero(i,2), SP_ECEF_ElevZero(i,3)] = geodetic2ecef(wgs84, SP_lla(i,1), SP_lla(i,2), 0.0);

    %% Plotting DEM
    
    if tag.plotDEM == 1
        %show DEM in lat-lon
        figDEM  = figure('Name', 'DEM', 'NumberTitle','off');
        mapshow(DEM_lla(:,:,2),DEM_lla(:,:,1),DEM_lla(:,:,3), 'DisplayType', 'surface')
        cbarDEM = colorbar;
        cbarDEM.Label.String = 'Surface elevation (m)';
        if DEMinfo.flatsurface==0
            caxis([min(min(DEM_lla(:,:,3))) max(max(DEM_lla(:,:,3)))])
        end
        colormap(parula), title('DEM in lat-lon')
        xlabel('Longitude (deg)'), ylabel('Latitude (deg)')
        saveas(figDEM, strjoin([nameOutputFile '_DEM.fig'],''));
%  if geoSYSp.Mode ~= "Single Point" , close (figDEM), end  % Plot only single point
    end

    %% -------- Compute the ionospheric delay at SP------------------------
    ionosphericDelay = ionoDelayComputation(ionisationLevelParam, GNSS_Week(i), GNSS_SoW(i),...
                                                RxECEF(i,:),TxECEF(i,:),SP_lla(i,:),TRAp);

    %% -------> Converting all coordinates from ECEF to ENU <--------------
    %  Note: the origin of ENU is in the Specular Point projected on the ellipsoid
    %  The rotation of ENU is enforced to have RX and TX in the XZ plane,
    %  forcing the Tx to be on the negative brunch of the axis
    %
    % -------   WGS - 84  (Local Normal direction = Z_axis) ---------------
    % -------   Geodetic Zenith (Normal to Geoidal Horizon) ---------------
    % -------   ENU is the local vertical reference system  ---------------
    % ------->>>Local reference System is rotated wrt ENU<<<--------------- 


    %Definign the current origin of the ENU reference system on the Specular Point
    %(the output is a col vector)
    %OENUinECEF = SP_ECEF_ElevZero(i,:).';
    OENUinECEF = SP_ECEF_meanElev(i,:).';
    %Including the OENU in TRAp
    TRAp.OENUinECEF = OENUinECEF;
    
    %% No error variables in TRAp (Amir)
    TRAp_NoError.OENUinECEF = OENUinECEF;

    TxENU  = ecef2enu_dc(TxECEF(i,:), OENUinECEF); 
    RxENU  = ecef2enu_dc(RxECEF(i,:), OENUinECEF);

    %The speed vector must not be traslated - just add an additional input
    VTxENU = ecef2enu_dc(VTxECEF(i,:), OENUinECEF, 'v'); 
    VRxENU = ecef2enu_dc(VRxECEF(i,:), OENUinECEF, 'v'); 
    
    %Checking that SP in ENU is zero - as it must be by definition
    SPENU = ecef2enu_dc(OENUinECEF, OENUinECEF);    
    if SPENU == [0,0,0], else, disp('ENU origin is wrong'), return, end  
        
    %azimuth angle of Tx in ENU
    phiTx = atan2(TxENU(2), TxENU(1));   
    %azimuth angle of Rx in ENU
    phiRX = atan2(RxENU(2), RxENU(1));     
    
    %%
    
    phi_Align_error = pi - abs(phiTx) - abs(phiRX);
    
    if phi_Align_error >= 0.01   % WARNING MESSAGE for MISALIGNMENT !
        Warning = strcat('Misalign(degs)=',num2str(phi_Align_error*180/pi));
        warndlg('GPS---SpecularPoint---Receiver',Warning)
    end        
    
    %ENU2Local: Rotation angle over the xy plane (i.e., azimuth) in RADians
    %pi- phiGPS is enforced to have the Tx always over the negative branch
    %of the x axis of the local reference system
    phi_ENU2local = pi - phiTx;
    
    TRAp.phi_ENU2local = phi_ENU2local;

    %% No error variables in TRAp (Amir)
    TRAp_NoError.phi_ENU2local = phi_ENU2local;

    %% AMIR TEST Rx Direction
    vRxLOCAL = ecef2local_dc_ak(TRAp.phi_ENU2local, VRxECEF(i,:), OENUinECEF,'v');
    Rx_heading_Local_rad = atan2(vRxLOCAL(2), vRxLOCAL(1));
    TRAp.Rx_heading_Local_deg = rad2deg(Rx_heading_Local_rad);

    vTxLOCAL = ecef2local_dc_ak(TRAp.phi_ENU2local, VTxECEF(i,:), OENUinECEF,'v');
    Tx_heading_Local_rad = atan2(vTxLOCAL(2), vTxLOCAL(1));
    TRAp.Tx_heading_Local_deg = rad2deg(Tx_heading_Local_rad);

     %% No error variables in TRAp (Amir)
    TRAp_NoError.Tx_heading_Local_deg = TRAp.Tx_heading_Local_deg;
    TRAp_NoError.Rx_heading_Local_deg = TRAp.Rx_heading_Local_deg;

    %% END AMIR TEST Rx Direction
 
    %-----ECEF TRA Parameters - Enforced as col vector -----%
    TRAp.TxECEF  = TxECEF(i,:).';
    TRAp.RxECEF  = RxECEF(i,:).';
    TRAp.VTxECEF = VTxECEF(i,:).';
    TRAp.VRxECEF = VRxECEF(i,:).';


    TRAp.Rx_roll = geoSYSp.Rx_roll(:,i);
    TRAp.Rx_pitch = geoSYSp.Rx_pitch(:,i);
    TRAp.Rx_yaw = geoSYSp.Rx_yaw(:,i);

    %% No error variables in TRAp (Amir)
    TRAp_NoError.TxECEF  = TxECEF(i,:).';
    TRAp_NoError.RxECEF  = RxECEF(i,:).';
    TRAp_NoError.VTxECEF = VTxECEF(i,:).';
    TRAp_NoError.VRxECEF = VRxECEF(i,:).';

    TRAp_NoError.Rx_roll = geoSYSp.Rx_roll(:,i);
    TRAp_NoError.Rx_pitch = geoSYSp.Rx_pitch(:,i);
    TRAp_NoError.Rx_yaw = geoSYSp.Rx_yaw(:,i);


%     if geoSYSp.Mode == "Trackwise" || geoSYSp.Mode == "Single Point"
%         TRAp.Rx_roll = geoSYSp.Rx_roll(:,i);
%         TRAp.Rx_pitch = geoSYSp.Rx_pitch(:,i);
%         TRAp.Rx_yaw = geoSYSp.Rx_yaw(:,i);
%     elseif  geoSYSp.Mode == "Monte Carlo" || geoSYSp.Mode == "Global"
%         TRAp.Rx_roll = geoSYSp.Rx_roll(i,:);
%         TRAp.Rx_pitch = geoSYSp.Rx_pitch(i,:);
%         TRAp.Rx_yaw = geoSYSp.Rx_yaw(i,:);
%     end


    %% Convertion of Tx and RX coordinates in rotated ENU: from now on LOCAL reference system (also SGR: SAVERS Global Reference)    
    %  Note: the rotation is clockwise, i.e., downwards
    
    %%These are column vector - as needed when applying transformations
    TRAp.Txlocal  = rotation_clockwise(TxENU,   phi_ENU2local, 'rad');
    TRAp.Rxlocal  = rotation_clockwise(RxENU,   phi_ENU2local, 'rad');
    TRAp.VTxlocal = rotation_clockwise(VTxENU,  phi_ENU2local, 'rad');
    TRAp.VRxlocal = rotation_clockwise(VRxENU,  phi_ENU2local, 'rad');
    %TRAp.SPlocal = rotation_clockwise([0,0,0], phi_ENU2local, 'rad');  

    %% No error variables in TRAp (Amir)
    TRAp_NoError.Txlocal = TRAp.Txlocal;
    TRAp_NoError.Rxlocal = TRAp.Rxlocal;
    TRAp_NoError.VTxlocal = TRAp.VTxlocal;
    TRAp_NoError.VRxlocal = TRAp.VRxlocal;
    
    %% --------------COMPUTATION OF INCIDENCE ANGLE------------------
    %  -----measured wrt the z axis, i.e., is a zenithal angle-------    
    %  ----An ad-hoc power point has been generated to explain this--

    %elevation angle in radians
     %elevation angle in radians
    gamma = -atan(TRAp.Txlocal(3)/TRAp.Txlocal(1));

    %incidence angle of Tx signal in local reference system
    incAngle_enuRot          = 90 - rad2deg(gamma);
    incAngle_enuRot_vect(i)  = incAngle_enuRot;
    TRAp.incAngle_enuRot     = incAngle_enuRot;

    %% No error variables in TRAp (Amir)
    TRAp_NoError.incAngle_enuRot = TRAp.incAngle_enuRot;

    %elevation angle - Alternatives
    incAngle_enuRot2 = atan2d(TRAp.Txlocal(3),TRAp.Txlocal(1));
    incAngle_enuRot3 = incAngle_enuRot + 90;    
    if round(incAngle_enuRot2,4) == round(incAngle_enuRot3,4), else, disp('Incidence angle is wrong'), return, end

%% ----  DEFINE THE SIZE OF THE INTEGRATION AREA ---------------

    %major semiaxis of the ellips at delay Dly_RNGmax (LEIMON report pp.11)
    % units are: [m]
    %Note: in Monte Carlo/Coarse Mode we minimize this size

    if      tag.CoarseMode == 0
            geoSYSp.Xmax_Range =  round(sqrt(2*2.99792458e2*(geoSYSp.Dly_RNGmax)*abs(TRAp.Rxlocal(3))*cosd(TRAp.incAngle_enuRot))/cosd(TRAp.incAngle_enuRot)^2/1000)*1000;
    elseif  tag.CoarseMode == 1
            geoSYSp.Xmax_Range = geoSYSp.AsideMC/2;
    end    

    % units are: [m]
    geoSYSp.Xmin_Range   = -geoSYSp.Xmax_Range; 
    % units are: [m]
    geoSYSp.Yabs_Range   =  geoSYSp.Xmax_Range;  
    
    %% ------------CONVERSION OF DEM FROM ECEF TO ENU ROTATED OF phi_ENU2local---------
    %waitbarEnuRot = waitbar(0,'MAPS FROM LAT-LON TO SAVERS ROTATED ENU');

    %%%Matrix implementation by DC
    pointDEM_enu           = ecef2enu_dc_mat(coordDEM_ecef, OENUinECEF);
    dataDEM_enuRot(:,:,1)  = pointDEM_enu(:,:,1).*cos(phi_ENU2local) - pointDEM_enu(:,:,2).*sin(phi_ENU2local); %coord x in rotated ENU of DEM point
    dataDEM_enuRot(:,:,2)  = pointDEM_enu(:,:,1).*sin(phi_ENU2local) + pointDEM_enu(:,:,2).*cos(phi_ENU2local); %coord y in rotated ENU of DEM point
    dataDEM_enuRot(:,:,3)  = pointDEM_enu(:,:,3);     

    %avoid the Earth curvature in case of simulation with active flat surface flag
    if DEMinfo.flatsurface==1
        dataDEM_enuRot(:,:,3)=0.0;
    end

    clear coordDEM_ecef

%     if tag.plotDEM == 1
%         %show DEM over the simulation area
%         figDEMsimArea = figure('Name', 'DEM in rotated ENU over the simulation area', 'NumberTitle','off');
%         mapshow(dataDEM_enuRot(:,:,1),dataDEM_enuRot(:,:,2),dataDEM_enuRot(:,:,3), 'DisplayType', 'surface')
%         axis([geoSYSp.Xmin_Range geoSYSp.Xmax_Range -geoSYSp.Yabs_Range  geoSYSp.Yabs_Range])        
%         cbarDEM              = colorbar;
%         cbarDEM.Label.String = 'Surface elevation (m)';
%         if DEMinfo.flatsurface==0
%             caxis([min(min(dataDEM_enuRot(:,:,3))) max(max(dataDEM_enuRot(:,:,3)))])
%         end
%         colormap(parula), title('DEM in rotated ENU over the simulation area')
%         xlabel('X axis (m)'), ylabel('Y axis (m)')
% 
%         if tag.savefigure == 1
%             saveas(figDEMsimArea, [dest_dir char(namefile) '_DEMinENUROTsimArea.fig']);
%         end
%     end
        
    %% -------> OPEN AND PROCESS THE LAND COVER MAP <--------------
    if size(coverMap_IDs,1) > 1

        %% CCI land cover map is referred to WGS84
        %elevationLandCoverMap = zeros(size(landCoverMap(:,:,1)));
        %%convert the cover map coordinates
        %[landCoverMap_ecef(:,:,1), landCoverMap_ecef(:,:,2), landCoverMap_ecef(:,:,3)] = geodetic2ecef(deg2rad(double(landCoverMap(:,:,1))),deg2rad(double(landCoverMap(:,:,2))),elevationLandCoverMap,wgs84);
        %landCoverMap_enu           = ecef2enu_dc_mat(landCoverMap_ecef, OENUinECEF);
        %coverMap_enuRot(:,:,1)  = landCoverMap_enu(:,:,1).*cos(phi_ENU2local) - landCoverMap_enu(:,:,2).*sin(phi_ENU2local); %coord x in rotated ENU of DEM point
        %coverMap_enuRot(:,:,2)  = landCoverMap_enu(:,:,1).*sin(phi_ENU2local) + landCoverMap_enu(:,:,2).*cos(phi_ENU2local); %coord y in rotated ENU of DEM point
        %coverMap_enuRot(:,:,3)  = landCoverMap(:,:,3);

        interpCoverMap = interp2(landCoverMap(:,:,2), landCoverMap(:,:,1), landCoverMap(:,:,3), DEM_lla(:,:,2), DEM_lla(:,:,1), 'nearest',-1);

        %convertion of the cover map from lat-lon to rotated ENU
        coverMap_enuRot        = dataDEM_enuRot;
        coverMap_enuRot(:,:,3) = interpCoverMap;
        coverMap_enuRot_all        = dataDEM_enuRot; %added by amir
        coverMap_enuRot_all(:,:,3) = interpCoverMap; %added by amir
    else
        coverMap_enuRot        = dataDEM_enuRot;
        coverMap_enuRot(:,:,3) = coverMap_IDs;

        interpCoverMap_all = interp2(landCoverMap_all(:,:,2), landCoverMap_all(:,:,1), landCoverMap_all(:,:,3), DEM_lla(:,:,2), DEM_lla(:,:,1), 'nearest',-1);

        %convertion of the cover map from lat-lon to rotated ENU
        coverMap_enuRot_all        = dataDEM_enuRot; %added by amir
        coverMap_enuRot_all(:,:,3) = interpCoverMap_all; %added by amir
    end

    clear DEM_lla
    %close(waitbarEnuRot) %======close waitbar

    %% ------------ TX ANTENNA GAIN AT SP -------------------------------------

    [TRAp.Yaw_TX, TRAp.TxGainAtSP,TRAp.TxAttitude_angles(3)]  = TxGain_SP(TRAp, GNSS_Week(i), GNSS_SoW(i));
    TRAp.TxGainAtSP_vect(i,:)       = TRAp.TxGainAtSP;
    TRAp.Yaw_TX_error = TRAp.TxAttitude_angles(3);

    %% No error variables in TRAp (Amir)
    [TRAp_NoError.Yaw_TX, TRAp_NoError.TxGainAtSP,TRAp_NoError.TxAttitude_angles(3)]  = TxGain_SP(TRAp_NoError, GNSS_Week(i), GNSS_SoW(i));
    TRAp_NoError.TxGainAtSP_vect(i,:)       = TRAp_NoError.TxGainAtSP;

    %% No error variables in geoSYSp (Amir)
    geoSYSp_NoError.Xmax_Range = geoSYSp.Xmax_Range;
    geoSYSp_NoError.Xmin_Range = geoSYSp.Xmin_Range;
    geoSYSp_NoError.Yabs_Range = geoSYSp.Yabs_Range;


    %% ----SIMULATION REGULAR GRID, ANGLES IN FACET REFERENCE SYSTEM, DELAY AND DOPPLER--------------

    [Xplot,Yplot, regGridSim_coverMap, sigmaParam, DDMparam, gammaParam] = computeGridCell_dc(TRAp, dataDEM_enuRot, coverMap_enuRot, ...
                                                                                 geoSYSp, nameOutputFile, tag, DEMinfo, ionosphericDelay,coverMap_enuRot_all);    
    

    %% No error (Amir)
    [~,~, ~, sigmaParam_NoError, DDMparam_NoError, gammaParam_NoError] = computeGridCell_dc(TRAp_NoError, dataDEM_enuRot, coverMap_enuRot, ...
                                                                                 geoSYSp_NoError, nameOutputFile, tag, DEMinfo, ionosphericDelay,coverMap_enuRot_all);    
    
%     path_factor_sp = gammaParam.PATH_factor_SP;
%     outputSimulations.PATH_factor(i,:) = gammaParam.PATH_factor_SP; %AMIR




    %number of grid cells in X
    Nx_surf = length(Xplot);
    %number of grid cells in Y
    Ny_surf = length(Yplot);

    surfParam.Nx_surf = Nx_surf; 
    surfParam.Ny_surf = Ny_surf; 
    surfParam.Xplot   = Xplot; 
    surfParam.Yplot   = Yplot; 
    
    %% ---- check cover classes present in the simulation area and       ---
    %  ---- set soil and vegetation parameters for each land cover class ---

    %number of classes in the integration area
    coverClassSimArea_IDs  = unique(regGridSim_coverMap);
    totNumOfClassesSimArea = size(coverClassSimArea_IDs,1);
%     if tag.directsignalonly == 0
    for l = 1:totNumOfClassesSimArea
        indexClass                       = find(coverMap_IDs == coverClassSimArea_IDs(l));
        terrainParameters_enuRot(l,:)    = terrainParameters(indexClass,:);
        vegetationParameters_enuRot(l,:) = vegetationParameters(indexClass,:);
    end
%     else
%     terrainParameters_enuRot = [0 1.0000 5.0000 3.0000 15.0000 1.1800 20.0000];
%     vegetationParameters_enuRot =[0 0 0 0 0 0 0 0 0 0 0 0];

    %Non ricordo più perchè qui si cerca, il sistema di rif è centrato su %SP
    indx_origin = Xplot == 0.0;
    indy_origin = Yplot == 0.0;

    TRAp.indx_origin = indx_origin;
    TRAp.indy_origin = indy_origin;

    %% No error variables in TRAp (Amir)
    TRAp_NoError.indx_origin = TRAp.indx_origin;
    TRAp_NoError.indy_origin = TRAp.indy_origin;

    %cover class at origin, at SP
%     if tag.directsignalonly == 0
    infoAtSPforRefFile.coverClassAtSP(i) = regGridSim_coverMap(indy_origin, indx_origin);
%     else
%     infoAtSPforRefFile.coverClassAtSP(i) = 3;   
%     end
    for indexTxFrequency = 1:length(TRAp.TxFrequency)
    %% --------- SIGMA ZERO COMPUTATION -----------------------------------
    % the sigma0 computation is carried out for each cover class, then the map 
    % of sigma0 over the simulation area is obtained assigning to each facet the sigma0
    % computed for the cover class at that facet.

        sigmaRR_cohe = zeros(Ny_surf,Nx_surf);
        sigmaRR_inco = zeros(Ny_surf,Nx_surf);
        sigmaLR_cohe = zeros(Ny_surf,Nx_surf);
        sigmaLR_inco = zeros(Ny_surf,Nx_surf);

        sigmaRR_cohe_NoError = zeros(Ny_surf,Nx_surf);
        sigmaRR_inco_NoError = zeros(Ny_surf,Nx_surf);
        sigmaLR_cohe_NoError = zeros(Ny_surf,Nx_surf);
        sigmaLR_inco_NoError = zeros(Ny_surf,Nx_surf);
 
        for m = 1:totNumOfClassesSimArea

            % check for unfeasable condition in forest: biomass zero and LAI
            % larger than zero
            if vegetationParameters_enuRot(m,2) == 1 && vegetationParameters_enuRot(m,7)==0.0 && vegetationParameters_enuRot(m,3)>0.0
                vegetationParameters_enuRot(m,3) = 0.0;
                disp('Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.')
                disp('LAI will be forced to zero.')
            end

            %these steps are used only when the sigma0 of the bare classes are
            %computed. These steps are necessary to avoid the computation of sigma0 over
            %the classes that are not bare 
            if vegetationParameters_enuRot(m,1) == 0 || vegetationParameters_enuRot(m,3)==0 && vegetationParameters_enuRot(m,7)==0
                sigmaParam.TH_s_flrIF = sigmaParam.TH_s_flr;
                sigmaParam.TH_s_flrIF(regGridSim_coverMap ~= coverClassSimArea_IDs(m)) = 9999;

                sigmaParam_NoError.TH_s_flrIF = sigmaParam_NoError.TH_s_flr;
                sigmaParam_NoError.TH_s_flrIF(regGridSim_coverMap ~= coverClassSimArea_IDs(m)) = 9999;
            end

            [sigmaRR_cohe_covIDm,sigmaLR_cohe_covIDm,sigmaRR_inco_covIDm,sigmaLR_inco_covIDm] =...
                    callFRoutine_IHS(Nx_surf, Ny_surf,TRAp.TxBeamwidth, TRAp.TxFrequency(indexTxFrequency),...
                                     TRAp.RxBeamwidth, TRAp.incAngle_enuRot, sigmaParam,...
                                     terrainParameters_enuRot(m,:), vegetationParameters_enuRot(m,:),TRAp,gammaParam,surfParam,geoSYSp);    
            if totNumOfClassesSimArea == 1
                sigmaRR_cohe = sigmaRR_cohe_covIDm;
                sigmaRR_inco = sigmaRR_inco_covIDm;
                sigmaLR_cohe = sigmaLR_cohe_covIDm;
                sigmaLR_inco = sigmaLR_inco_covIDm;
            else
                sigmaRR_cohe(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaRR_cohe_covIDm(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaRR_inco(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaRR_inco_covIDm(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaLR_cohe(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaLR_cohe_covIDm(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaLR_inco(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaLR_inco_covIDm(regGridSim_coverMap == coverClassSimArea_IDs(m));
            end

                    [sigmaRR_cohe_covIDm_NoError,sigmaLR_cohe_covIDm_NoError,sigmaRR_inco_covIDm_NoError,sigmaLR_inco_covIDm_NoError] =...
                    callFRoutine_IHS(Nx_surf, Ny_surf,TRAp_NoError.TxBeamwidth, TRAp_NoError.TxFrequency(indexTxFrequency),...
                                     TRAp_NoError.RxBeamwidth, TRAp_NoError.incAngle_enuRot, sigmaParam_NoError,...
                                     terrainParameters_enuRot(m,:), vegetationParameters_enuRot(m,:),TRAp_NoError,gammaParam_NoError,surfParam,geoSYSp_NoError);    
            if totNumOfClassesSimArea == 1
                sigmaRR_cohe_NoError = sigmaRR_cohe_covIDm_NoError;
                sigmaRR_inco_NoError = sigmaRR_inco_covIDm_NoError;
                sigmaLR_cohe_NoError = sigmaLR_cohe_covIDm_NoError;
                sigmaLR_inco_NoError = sigmaLR_inco_covIDm_NoError;
            else
                sigmaRR_cohe_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaRR_cohe_covIDm_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaRR_inco_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaRR_inco_covIDm_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaLR_cohe_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaLR_cohe_covIDm_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m));
                sigmaLR_inco_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m)) = sigmaLR_inco_covIDm_NoError(regGridSim_coverMap == coverClassSimArea_IDs(m));
            end           


        end
        
        %% ---- plot the scattering diagrams ---
        if tag.plotSigma0 == 1

            sigmaLR_inco_plot = sigmaLR_inco;
            sigmaLR_cohe_plot = sigmaLR_cohe;
            sigmaRR_inco_plot = sigmaRR_inco;
            sigmaRR_cohe_plot = sigmaRR_cohe;


            figScatDiagTot=figure('Name', 'Bistatic scattering coefficient diagram', 'NumberTitle','off');
            set(gcf,'Position',[0 50 900 700])

            subplot(2,2,3) % lower left
            contourf(Xplot,Yplot,sigmaLR_inco_plot(:,:,1),'LineColor', 'none');
            cbarDiag1=colorbar;
            cbarDiag1.Label.String = '\sigma_0 (m^2/m^2)';
            if indexTxFrequency==1
            title(' LR1 incoherent ')
            elseif indexTxFrequency==2
            title(' LR5 incoherent ')
            end
            xlabel('X axis (m)')
            ylabel('Y axis (m)')

            subplot(2,2,1) % upper left
            contourf(Xplot,Yplot,sigmaLR_cohe_plot(:,:,1),'LineColor', 'none');
            cbarDiag1=colorbar;
            cbarDiag1.Label.String = '\sigma_0 (m^2/m^2)';
            if indexTxFrequency==1
            title(' LR1 coherent ')
            elseif indexTxFrequency==2
            title(' LR5 coherent ')
            end
            xlabel('X axis (m)')
            ylabel('Y axis (m)')

            subplot(2,2,4) % lower right
            contourf(Xplot,Yplot,sigmaRR_inco_plot(:,:,1),'LineColor', 'none');
            cbarDiag1=colorbar;
            cbarDiag1.Label.String = '\sigma_0 (m^2/m^2)';
            if indexTxFrequency==1
            title(' RR1 incoherent ')
            elseif indexTxFrequency==2
            title(' RR5 incoherent ')
            end
            xlabel('X axis (m)')
            ylabel('Y axis (m)')

            subplot(2,2,2) % upper right
            contourf(Xplot,Yplot,sigmaRR_cohe_plot(:,:,1),'LineColor', 'none');
            cbarDiag1=colorbar;
            cbarDiag1.Label.String = '\sigma_0 (m^2/m^2)';
            if indexTxFrequency==1
            title(' RR1 coherent ')
            elseif indexTxFrequency==2
            title(' RR5 coherent ')
            end
            xlabel('X axis (m)')
            ylabel('Y axis (m)')
            
            saveas(figScatDiagTot, strjoin([nameOutputFile '_scatDiag.fig'],''));
        end          
                    
        %% --------- COMPUTE DDM AND PEAK REFLECTIVITY ---------------------------
        %this routine integrates the scattering of all coverages
        %and compute DDM and maximum of DDM
        
        %% Set random seed for all random generator during process for Noise (for MM)
        randseed = rand;
        geoSYSp.randseed = floor(randseed * 1e6);
        geoSYSp_NoError.randseed = geoSYSp.randseed;

        [max_sigPower_noNoise_LR(i,indexTxFrequency), max_sigPower_noNoise_RR(i,indexTxFrequency), max_sigPower_Noise_LR(i,indexTxFrequency), max_sigPower_Noise_RR(i,indexTxFrequency),...
                sigPower_noNoise_LR_lin(:,:,i,indexTxFrequency), sigPower_noNoise_RR_lin(:,:,i,indexTxFrequency), sigPower_Noise_LR_lin(:,:,i,indexTxFrequency), sigPower_Noise_RR_lin(:,:,i,indexTxFrequency),...
                NLbb_DOWN_LR(:,:,i,indexTxFrequency),NLbb_DOWN_RR(:,:,i,indexTxFrequency), LR_UP_Power(i,indexTxFrequency), RR_UP_Power(i,indexTxFrequency), LR_UP_Power_Flu(i,indexTxFrequency), RR_UP_Power_Flu(i,indexTxFrequency),...
                PN_star(i,indexTxFrequency), NLbb_UP(:,:,i,indexTxFrequency), UP_pattern(i,indexTxFrequency), TxGainAtRx(i,indexTxFrequency), geoSYSp, AcohLR(i,indexTxFrequency), AdiffLR(i,indexTxFrequency), rho_TRsp(i,indexTxFrequency)] =...
                    BISTATIC_IHS_dc(Ny_surf, Nx_surf, sigmaRR_cohe, sigmaRR_inco, sigmaLR_cohe, sigmaLR_inco,...
                                    DDMparam, geoSYSp, TRAp, nameOutputFile, tag, surfParam, instrumentConfParam, indexTxFrequency);
         

        [max_sigPower_noNoise_NoError_LR(i,indexTxFrequency), max_sigPower_noNoise_NoError_RR(i,indexTxFrequency), max_sigPower_Noise_NoError_LR(i,indexTxFrequency), max_sigPower_Noise_NoError_RR(i,indexTxFrequency),...
                sigPower_noNoise_NoError_LR_lin(:,:,i,indexTxFrequency), sigPower_noNoise_NoError_RR_lin(:,:,i,indexTxFrequency), sigPower_Noise_NoError_LR_lin(:,:,i,indexTxFrequency), sigPower_Noise_NoError_RR_lin(:,:,i,indexTxFrequency),...
                NLbb_DOWN_NoError_LR(:,:,i,indexTxFrequency),NLbb_DOWN_NoError_RR(:,:,i,indexTxFrequency), LR_UP_Power_NoError(i,indexTxFrequency), RR_UP_Power_NoError(i,indexTxFrequency), LR_UP_Power_Flu_NoError(i,indexTxFrequency), RR_UP_Power_Flu_NoError(i,indexTxFrequency),...
                PN_star_NoError(i,indexTxFrequency), NLbb_UP_NoError(:,:,i,indexTxFrequency), UP_pattern_NoError(i,indexTxFrequency), TxGainAtRx_NoError(i,indexTxFrequency), geoSYSp_NoError, AcohLR_NoError(i,indexTxFrequency), AdiffLR_NoError(i,indexTxFrequency), rho_TRsp_NoError(i,indexTxFrequency)] =...
                    BISTATIC_IHS_dc(Ny_surf, Nx_surf, sigmaRR_cohe_NoError, sigmaRR_inco_NoError, sigmaLR_cohe_NoError, sigmaLR_inco_NoError,...
                                    DDMparam_NoError, geoSYSp_NoError, TRAp_NoError, nameOutputFile, tag, surfParam, instrumentConfParam_NoError, indexTxFrequency);

        %if trackwise mode then close the figures
        if Nruns>1
            close all
        end

%         %% AMIR
%         ponsurf_temp = reflcomponents(i).PonSurf;
%         ponsurf_temp(isnan(ponsurf_temp)) = 0;
%         ponsurf = db(max(ponsurf_temp(:)), 'power');  % Flatten with (:) for max
%         
%         centre_x = ceil(Nx_surf / 2);
%         centre_y = ceil(Ny_surf / 2);
%         gain_temp = reflcomponents(i).Gain_LHCP_onsurf_lin(centre_x, centre_y);
%         PATH_factor_temp = reflcomponents(i).PATH_factor(centre_x, centre_y);
%         
%         outputSimulations.reflcomponentx.max_sigPower_noNoise_LR(i, :) = reflcomponents(i).max_sigPower_noNoise_LR;
%         outputSimulations.reflcomponentx.max_sigPower_Noise_LR(i, :)   = reflcomponents(i).max_sigPower_Noise_LR;
%         outputSimulations.reflcomponentx.TxGainAtRx(i, :)              = reflcomponents(i).TxGainAtRx;
%         outputSimulations.reflcomponentx.PonSurf(i, :)                 = ponsurf;
%         outputSimulations.reflcomponentx.LR_UP_Power(i, :)             = reflcomponents(i).LR_UP_Power;
%         outputSimulations.reflcomponentx.Gain_LHCP_onsurf_lin(i, :)    = gain_temp;
%         outputSimulations.reflcomponentx.PATH_factor(i, :)             = PATH_factor_temp;
% % Individual components of reflectivity_noNoise_LR_tot
% component_sigPower_noNoise   = sigPower_noNoise_LR_lin(:, :, i, indexTxFrequency);
% component_4pi_squared        = (4 * pi)^2;
% component_squaredDistanceSum = gammaParam.squaredDistanceSum;
% component_lambda_squared     = TRAp.waveL(indexTxFrequency)^2;
% component_Gain_LHCP          = gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency);
% component_TxSignalPower      = db2pow(TRAp.TxSignalPower(indexTxFrequency));
% component_TxGainAtSP         = TRAp.TxGainAtSP(indexTxFrequency);
% 
% % Final reflectivity calculation
% reflectivity_noNoise_LR_tot = ...
%     component_sigPower_noNoise .* ...
%     component_4pi_squared .* ...
%     component_squaredDistanceSum ./ ...
%     (component_lambda_squared * component_Gain_LHCP * component_TxSignalPower * component_TxGainAtSP);
% 
% % Save all components into outputSimulations.reflcomponentx
% % outputSimulations.reflcomponentx.sigPower_noNoise_LR_lin(i, :)   = component_sigPower_noNoise;
% outputSimulations.reflcomponentx.fourPi_squared(i, :)            = component_4pi_squared;
% outputSimulations.reflcomponentx.squaredDistanceSum(i, :)        = component_squaredDistanceSum;
% outputSimulations.reflcomponentx.lambda_squared(i, :)            = component_lambda_squared;
% outputSimulations.reflcomponentx.Gain_LHCP_onsurf_atSP(i, :)     = component_Gain_LHCP;
% outputSimulations.reflcomponentx.TxSignalPower_lin(i, :)         = component_TxSignalPower;
% outputSimulations.reflcomponentx.TxGainAtSP(i, :)                = component_TxGainAtSP;
% refltemp_lin = max(max(sigPower_noNoise_LR_lin(:,:,i,indexTxFrequency)));
% outputSimulations.reflcomponentx.sigPower_noNoise_LR_lin(i,:)= refltemp_lin;
% % outputSimulations.reflcomponentx.reflectivity_noNoise_LR_tot(i, :)  = reflectivity_noNoise_LR_tot;
% 
% 
%         %% AMIR      %% ----------- REFLECTIVITY COMPUTATION -------------------------------

        %-----REFLECTIVITY FROM NOISE FREE POWER------

%         reflectivity_noNoise_LR_tot             = sigPower_noNoise_LR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
%         reflectivity_noNoise_LR_tot_dB          = db(reflectivity_noNoise_LR_tot,'power');
    
        reflectivity_noNoise_NoError_LR_tot             = sigPower_noNoise_NoError_LR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
        reflectivity_noNoise_NoError_LR_tot_dB          = db(reflectivity_noNoise_NoError_LR_tot,'power');

%         reflectivity_noNoise_RR_tot             = sigPower_noNoise_RR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP((indexTxFrequency)));
%         reflectivity_noNoise_RR_tot_dB          = db(reflectivity_noNoise_RR_tot,'power');
%         
        reflectivity_noNoise_NoError_RR_tot             = sigPower_noNoise_NoError_RR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP((indexTxFrequency)));
        reflectivity_noNoise_NoError_RR_tot_dB          = db(reflectivity_noNoise_NoError_RR_tot,'power');

        if tag.saveReflectivityMat==1 && indexTxFrequency == 1
            save(strjoin([nameOutputFile '_reflectivityLIN_LR.mat'],''), 'reflectivity_noNoise_NoError_LR_tot');
            save(strjoin([nameOutputFile '_reflectivityLIN_RR.mat'],''), 'reflectivity_noNoise_NoError_RR_tot'); 
        elseif tag.saveReflectivityMat==1 && indexTxFrequency == 2
            save(strjoin([nameOutputFile '_reflectivityLIN_LR_5.mat'],''), 'reflectivity_noNoise_NoError_LR_tot');
            save(strjoin([nameOutputFile '_reflectivityLIN_RR_5.mat'],''), 'reflectivity_noNoise_NoError_RR_tot'); 
        end

        reflectivity_Noisy_LR_tot             = sigPower_Noise_NoError_LR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
        reflectivity_Noisy_LR_tot_dB          = db(reflectivity_Noisy_LR_tot,'power');

        reflectivity_Noisy_RR_tot             = sigPower_Noise_NoError_RR_lin(:,:,i,indexTxFrequency).*(4*pi)^2.*gammaParam.squaredDistanceSum./(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP((indexTxFrequency)));
        reflectivity_Noisy_RR_tot_dB          = db(reflectivity_Noisy_RR_tot,'power');

%% Added by amir for GREAT - DDM of reflectivity

if tag.plotDDM == 1
Fplot = 1:1:geoSYSp.ny_Dfreq;
Rplot = 1:1:geoSYSp.nx_Rdelay;

    if isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
    figDDM_LRtotRefl = figure('Name', 'LR reflectivity clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,reflectivity_noNoise_NoError_LR_tot_dB);
    colorbar
    colormap('parula')
    axis xy
    caxis([max(reflectivity_noNoise_NoError_LR_tot_dB(:)) - 20, max(reflectivity_noNoise_NoError_LR_tot_dB(:))])

    if indexTxFrequency==1
    title(' LR1 reflectivity clean DDM in dB')
    elseif indexTxFrequency==2
    title(' LR5 reflectivity clean DDM in dB')
    end

    figDDM_LRtotRefl_lin = figure('Name', 'LR reflectivity clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(reflectivity_noNoise_NoError_LR_tot_dB));
    colorbar
    colormap('parula')
    axis xy

    if indexTxFrequency==1
    title(' LR1 reflectivity clean DDM')
    elseif indexTxFrequency==2
    title(' LR5 reflectivity clean DDM')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')
    caxis([0, max(db2pow(reflectivity_noNoise_NoError_LR_tot_dB(:)))])

     if indexTxFrequency==1
        saveas(figDDM_LRtotRefl, strjoin([nameOutputFile '_DDMnoiseFree_LR1_refl.fig'],''));
        saveas(figDDM_LRtotRefl_lin, strjoin([nameOutputFile '_DDMnoiseFree_LR1_refl_lin.fig'],''));
    else
        saveas(figDDM_LRtotRefl, strjoin([nameOutputFile '_DDMnoiseFree_LR5_refl.fig'],''));
        saveas(figDDM_LRtotRefl_lin, strjoin([nameOutputFile '_DDMnoiseFree_LR5_refl_lin.fig'],''));    
     end

    end


    if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    figDDM_LRtotReflNoisy = figure('Name', 'LR reflectivity noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,reflectivity_Noisy_LR_tot_dB);
    colorbar
    colormap('parula')
    axis xy
    caxis([max(reflectivity_Noisy_LR_tot_dB(:)) - 20, max(reflectivity_Noisy_LR_tot_dB(:))])

    if indexTxFrequency==1
    title(' LR1 reflectivity noisy DDM in dB')
    elseif indexTxFrequency==2
    title(' LR5 reflectivity noisy DDM in dB')
    end
    
    figDDM_LRtotReflNoisy_lin = figure('Name', 'LR reflectivity noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(reflectivity_Noisy_LR_tot_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(reflectivity_Noisy_LR_tot_dB(:)))])

    if indexTxFrequency==1
    title(' LR1 reflectivity noisy DDM')
    elseif indexTxFrequency==2
    title(' LR5 reflectivity noisy DDM')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')
    
    if indexTxFrequency==1
        saveas(figDDM_LRtotReflNoisy, strjoin([nameOutputFile '_DDMnoisy_LR1_refl.fig'],''));
        saveas(figDDM_LRtotReflNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_LR1_refl_lin.fig'],''));
    else
        saveas(figDDM_LRtotReflNoisy, strjoin([nameOutputFile '_DDMnoisy_LR5_refl.fig'],''));
        saveas(figDDM_LRtotReflNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_LR5_refl_lin.fig'],''));    
    end
    
    end


end

if tag.plotDDM == 1

    if isfield(TRAp, 'InternalFlag_NoError') && TRAp.InternalFlag_NoError == 1
    figDDM_RRtotRefl = figure('Name', 'RR reflectivity clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,reflectivity_noNoise_NoError_RR_tot_dB);
    colorbar
    colormap('parula')
    axis xy
    caxis([max(reflectivity_noNoise_NoError_RR_tot_dB(:)) - 20, max(reflectivity_noNoise_NoError_RR_tot_dB(:))])

    if indexTxFrequency==1
    title(' RR1 reflectivity clean DDM in dB')
    elseif indexTxFrequency==2
    title(' RR5 reflectivity clean DDM in dB')
    end

    figDDM_RRtotRefl_lin = figure('Name', 'RR reflectivity clean DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(reflectivity_noNoise_NoError_RR_tot_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(reflectivity_noNoise_NoError_RR_tot_dB(:)))])

    if indexTxFrequency==1
    title(' RR1 reflectivity clean DDM')
    elseif indexTxFrequency==2
    title(' RR5 reflectivity clean DDM')
    end

    xlabel('Delay bin')
    ylabel('Doppler bin')

    if indexTxFrequency==1
        saveas(figDDM_RRtotRefl, strjoin([nameOutputFile '_DDMnoiseFree_RR1_refl.fig'],''));
        saveas(figDDM_RRtotRefl_lin, strjoin([nameOutputFile '_DDMnoiseFree_RR1_refl_lin.fig'],''));
    else
        saveas(figDDM_RRtotRefl, strjoin([nameOutputFile '_DDMnoiseFree_RR5_refl.fig'],''));
        saveas(figDDM_RRtotRefl_lin, strjoin([nameOutputFile '_DDMnoiseFree_RR5_refl_lin.fig'],''));    
    end


    end
    
    if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
    figDDM_RRtotReflNoisy = figure('Name', 'RR reflectivity noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,reflectivity_Noisy_RR_tot_dB);
    colorbar
    colormap('parula')
    axis xy
    caxis([max(reflectivity_Noisy_RR_tot_dB(:)) - 20, max(reflectivity_Noisy_RR_tot_dB(:))])

    if indexTxFrequency==1
    title(' RR1 reflectivity noisy DDM in dB')
    elseif indexTxFrequency==2
    title(' RR5 reflectivity noisy DDM in dB')
    end

    figDDM_RRtotReflNoisy_lin = figure('Name', 'RR reflectivity noisy DDM', 'NumberTitle','off');
    imagesc(Rplot,Fplot,db2pow(reflectivity_Noisy_RR_tot_dB));
    colorbar
    colormap('parula')
    axis xy
    caxis([0, max(db2pow(reflectivity_Noisy_RR_tot_dB(:)))])

    xlabel('Delay bin')
    ylabel('Doppler bin')
    
    if indexTxFrequency==1
    title(' RR1 reflectivity noisy DDM')
    elseif indexTxFrequency==2
    title(' RR5 reflectivity noisy DDM')
    end

    if indexTxFrequency==1
        saveas(figDDM_RRtotReflNoisy, strjoin([nameOutputFile '_DDMnoisy_RR1_refl.fig'],''));
        saveas(figDDM_RRtotReflNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_RR1_refl_lin.fig'],''));
    else
        saveas(figDDM_RRtotReflNoisy, strjoin([nameOutputFile '_DDMnoisy_RR5_refl.fig'],''));
        saveas(figDDM_RRtotReflNoisy_lin, strjoin([nameOutputFile '_DDMnoisy_RR5_refl_lin.fig'],''));    
    end

    end
end

%% Added by amir for GREAT - DDM of reflectivity - End        
        %reflectivity at peak of DDM
%         max_reflectivity_noNoise_LR_tot_dB(i,indexTxFrequency)   = max(max(reflectivity_noNoise_LR_tot_dB));
%         max_reflectivity_noNoise_RR_tot_dB(i,indexTxFrequency)   = max(max(reflectivity_noNoise_RR_tot_dB));

        max_reflectivity_noNoise_NoError_LR_tot_dB(i,indexTxFrequency)   = max(max(reflectivity_noNoise_NoError_LR_tot_dB));
        max_reflectivity_noNoise_NoError_RR_tot_dB(i,indexTxFrequency)   = max(max(reflectivity_noNoise_NoError_RR_tot_dB));
        

        %reflectivity in a small window around the nominal SP
        
        geobWinWidthInDelay = geoSYSp.geoboundedWindowWidthInDelay_numBins;
        geobWinWidthInFreq  = geoSYSp.geoboundedWindowWidthInFreq_numBins;

%         geoboundedReflecDDM_noNoise_LR                   = reflectivity_noNoise_LR_tot_dB(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2));
%         max_geobReflec_noNoise_LR_dB(i,indexTxFrequency) = max(max(geoboundedReflecDDM_noNoise_LR));

        geoboundedReflecDDM_noNoise_NoError_LR                   = reflectivity_noNoise_NoError_LR_tot_dB(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2));
        max_geobReflec_noNoise_NoError_LR_dB(i,indexTxFrequency) = max(max(geoboundedReflecDDM_noNoise_NoError_LR));
        

%         geoboundedReflecDDM_noNoise_RR                   = reflectivity_noNoise_RR_tot_dB(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2));
%         max_geobReflec_noNoise_RR_dB(i,indexTxFrequency) = max(max(geoboundedReflecDDM_noNoise_RR));

        geoboundedReflecDDM_noNoise_NoError_RR                   = reflectivity_noNoise_NoError_RR_tot_dB(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2));
        max_geobReflec_noNoise_NoError_RR_dB(i,indexTxFrequency) = max(max(geoboundedReflecDDM_noNoise_NoError_RR));
        


        %-----REFLECTIVITY FROM NOISY POWER------

        %%estimate the noise in received signal, to compute the reflectivity of the
        %%noisy received signal

        % found a window of the DDM where there is no signal
        [~, columnIndexPixelsWithSignal]=find(sigPower_noNoise_LR_lin(:,:,i,indexTxFrequency));
        mincolumnPixelsWithSignal=min(columnIndexPixelsWithSignal);

        if mincolumnPixelsWithSignal==1
            %estimate the noise in the window of the noisy DDM where there is no signal
            noiseEstimateLR_bottomBox=mean(mean(sigPower_Noise_LR_lin(1:3,1:10,i,indexTxFrequency)));
            noiseEstimateLR_upBox=mean(mean(sigPower_Noise_LR_lin(geoSYSp.ny_Dfreq-3:geoSYSp.ny_Dfreq,1:5,i,indexTxFrequency)));
            noiseEstimateLR=(noiseEstimateLR_bottomBox+noiseEstimateLR_upBox)/2;
            noiseEstimateRR_bottomBox=mean(mean(sigPower_Noise_RR_lin(1:3,1:10,i,indexTxFrequency)));
            noiseEstimateRR_upBox=mean(mean(sigPower_Noise_RR_lin(geoSYSp.ny_Dfreq-3:geoSYSp.ny_Dfreq,1:5,i,indexTxFrequency)));
            noiseEstimateRR=(noiseEstimateRR_bottomBox+noiseEstimateRR_upBox)/2;
        else
            %estimate the noise in the window of the noisy DDM where there is no signal
            noiseEstimateLR=mean(mean(sigPower_Noise_LR_lin(:,(1:mincolumnPixelsWithSignal-1),i,indexTxFrequency)));
            noiseEstimateRR=mean(mean(sigPower_Noise_RR_lin(:,(1:mincolumnPixelsWithSignal-1),i,indexTxFrequency)));
        end

%             [~, columnIndexPixelsWithSignal_NoError]=find(sigPower_noNoise_NoError_LR_lin(:,:,i,indexTxFrequency));
%             mincolumnPixelsWithSignal_NoError=min(columnIndexPixelsWithSignal_NoError);
% 
%         if mincolumnPixelsWithSignal_NoError==1
%             %estimate the noise in the window of the noisy DDM where there is no signal
%             noiseEstimateLR_NoError_bottomBox=mean(mean(sigPower_Noise_NoError_LR_lin(1:3,1:10,i,indexTxFrequency)));
%             noiseEstimateLR_NoError_upBox=mean(mean(sigPower_Noise_NoError_LR_lin(geoSYSp.ny_Dfreq-3:geoSYSp.ny_Dfreq,1:5,i,indexTxFrequency)));
%             noiseEstimateLR_NoError=(noiseEstimateLR_NoError_bottomBox+noiseEstimateLR_NoError_upBox)/2;
%             noiseEstimateRR_NoError_bottomBox=mean(mean(sigPower_Noise_NoError_RR_lin(1:3,1:10,i,indexTxFrequency)));
%             noiseEstimateRR_NoError_upBox=mean(mean(sigPower_Noise_NoError_RR_lin(geoSYSp.ny_Dfreq-3:geoSYSp.ny_Dfreq,1:5,i,indexTxFrequency)));
%             noiseEstimateRR_NoError=(noiseEstimateRR_NoError_bottomBox+noiseEstimateRR_NoError_upBox)/2;
%         else
%             %estimate the noise in the window of the noisy DDM where there is no signal
%             noiseEstimateLR_NoError=mean(mean(sigPower_Noise_NoError_LR_lin(:,(1:mincolumnPixelsWithSignal_NoError-1),i,indexTxFrequency)));
%             noiseEstimateRR_NoError=mean(mean(sigPower_Noise_NoError_RR_lin(:,(1:mincolumnPixelsWithSignal_NoError-1),i,indexTxFrequency)));
%         end

        %remove the estimated noise from the peak of the noisy DDM    
        max_reflectivity_Noisy_LR_tot                        = ((max_sigPower_Noise_LR(i,indexTxFrequency)-noiseEstimateLR)/geoSYSp.G_sys_lin_LHCP_nominal(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
        max_reflectivity_Noisy_LR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_LR_tot,'power');
        max_reflectivity_Noisy_RR_tot                        = ((max_sigPower_Noise_RR(i,indexTxFrequency)-noiseEstimateRR)/geoSYSp.G_sys_lin_RHCP_nominal(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
        max_reflectivity_Noisy_RR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_RR_tot,'power');


%% Noisy Cal - the same approach as MM calibration by rongowai - By Amir
% 
% LL = gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency);
% RR = gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency);
% RL = gammaParam.Gain_RL_onsurf_lin_atSP(indexTxFrequency);
% LR = gammaParam.Gain_LR_onsurf_lin_atSP(indexTxFrequency);
% 
% detG = LL.*RR - RL.*LR;
% detG(abs(detG)<eps) = eps;
% 
% inv11 = RR./detG;
% inv12 = -RL./detG;
% inv21 = -LR./detG;
% inv22 = LL./detG;
% 
% P_L = double(max_sigPower_Noise_LR(i,indexTxFrequency) - noiseEstimateLR);
% P_R = double(max_sigPower_Noise_RR(i,indexTxFrequency) - noiseEstimateRR);
% 
% t3_L = ((4*pi)^2 * gammaParam.squaredDistanceSum) ./ ...
%         (TRAp.waveL(indexTxFrequency)^2 * ...
%          db2pow(TRAp.TxSignalPower(indexTxFrequency)) * ...
%          TRAp.TxGainAtSP(indexTxFrequency) * ...
%          geoSYSp.G_sys_lin_LHCP_nominal(indexTxFrequency));
% 
% t3_R = ((4*pi)^2 * gammaParam.squaredDistanceSum) ./ ...
%         (TRAp.waveL(indexTxFrequency)^2 * ...
%          db2pow(TRAp.TxSignalPower(indexTxFrequency)) * ...
%          TRAp.TxGainAtSP(indexTxFrequency) * ...
%          geoSYSp.G_sys_lin_RHCP_nominal(indexTxFrequency));
% 
% Rcop = t3_L .* (inv11 .* P_L + inv12 .* P_R);
% Rxpl = t3_R .* (inv21 .* P_L + inv22 .* P_R);
% 
% max_reflectivity_Noisy_LR_tot = max(Rcop);
% max_reflectivity_Noisy_RR_tot = max(Rxpl);
% 
% max_reflectivity_Noisy_LR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_LR_tot,'power');
% max_reflectivity_Noisy_RR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_RR_tot,'power');


%% END


%         max_reflectivity_Noisy_NoError_LR_tot                        = ((max_sigPower_Noise_NoError_LR(i,indexTxFrequency)-noiseEstimateLR_NoError)/geoSYSp.G_sys_lin_LHCP_nominal(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
%         max_reflectivity_Noisy_NoError_LR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_NoError_LR_tot,'power');
%         max_reflectivity_Noisy_NoError_RR_tot                        = ((max_sigPower_Noise_NoError_RR(i,indexTxFrequency)-noiseEstimateRR_NoError)/geoSYSp.G_sys_lin_RHCP_nominal(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
%         max_reflectivity_Noisy_NoError_RR_tot_dB(i,indexTxFrequency) = db(max_reflectivity_Noisy_NoError_RR_tot,'power');
%         

        %remove the estimated noise from the geobounded peak of the noisy DDM
        geoboundedNoisyPowerDDM_LR = sigPower_Noise_LR_lin(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2),i,indexTxFrequency);
        max_geobPower_Noisy_LR(i,indexTxFrequency) = max(max(geoboundedNoisyPowerDDM_LR));

%         geoboundedNoisyPowerDDM_NoError_LR = sigPower_Noise_NoError_LR_lin(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2),i,indexTxFrequency);
%         max_geobPower_Noisy_NoError_LR(i,indexTxFrequency) = max(max(geoboundedNoisyPowerDDM_NoError_LR));
%         
        %NOTE: in some cases the noise estimated over the full DDM is
        %higher than the max of the geobouded DDM, therefore it is not
        %possible to remove the noise frm the geobounded DDM
%         max_geobReflec_Noisy_LR    = ((max_geobPower_Noisy_LR-noiseEstimateLR)/geoSYSp.G_sys_lin_LHCP(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
%         max_geobReflec_Noisy_LR_dB(i,indexTxFrequency) = db(max_geobReflec_Noisy_LR,'power');

        geoboundedNoisyPowerDDM_RR = sigPower_Noise_RR_lin(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2),i,indexTxFrequency);
        max_geobPower_Noisy_RR(i,indexTxFrequency) = max(max(geoboundedNoisyPowerDDM_RR));

%         geoboundedNoisyPowerDDM_NoError_RR = sigPower_Noise_NoError_RR_lin(floor((geoSYSp.ny_Dfreq/2+1)-(geobWinWidthInFreq-1)/2):ceil((geoSYSp.ny_Dfreq/2+1)+(geobWinWidthInFreq-1)/2),floor((geoSYSp.delayOffset_bin+1)-(geobWinWidthInDelay(indexTxFrequency)-1)/2):ceil((geoSYSp.delayOffset_bin+1)+(geobWinWidthInDelay(indexTxFrequency)-1)/2),i,indexTxFrequency);
%         max_geobPower_Noisy_NoError_RR(i,indexTxFrequency) = max(max(geoboundedNoisyPowerDDM_NoError_RR));
%         
        %NOTE: in some cases the noise estimated over the full DDM is
        %higher than the max of the geobouded DDM, therefore it is not
        %possible to remove the noise frm the geobounded DDM
%         max_geobReflec_Noisy_RR    = ((max_geobPower_Noisy_RR-noiseEstimateRR)/geoSYSp.G_sys_lin_RHCP(indexTxFrequency))*(4*pi)^2*gammaParam.squaredDistanceSum/(TRAp.waveL(indexTxFrequency)^2*gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency)*db2pow(TRAp.TxSignalPower(indexTxFrequency))*TRAp.TxGainAtSP(indexTxFrequency));
%         max_geobReflec_Noisy_RR_dB(i,indexTxFrequency) = db(max_geobReflec_Noisy_RR,'power');

        %compute the SNR in dB
        if tag.SNR==1
            SNR_LR(i,indexTxFrequency) = db((max_sigPower_Noise_LR(i,indexTxFrequency)-noiseEstimateLR)/noiseEstimateLR,'power');
            SNR_RR(i,indexTxFrequency) = db((max_sigPower_Noise_RR(i,indexTxFrequency)-noiseEstimateRR)/noiseEstimateRR,'power');

%             SNR_NoError_LR(i,indexTxFrequency) = db((max_sigPower_Noise_NoError_LR(i,indexTxFrequency)-noiseEstimateLR_NoError)/noiseEstimateLR_NoError,'power');
%             SNR_NoError_RR(i,indexTxFrequency) = db((max_sigPower_Noise_NoError_RR(i,indexTxFrequency)-noiseEstimateRR_NoError)/noiseEstimateRR_NoError,'power');
%             
%             geobSNR_LR(i,indexTxFrequency) = db((max_geobPower_Noisy_LR-noiseEstimateLR)/noiseEstimateLR,'power');
%             geobSNR_RR(i,indexTxFrequency) = db((max_geobPower_Noisy_RR-noiseEstimateRR)/noiseEstimateRR,'power');

        end


        infoAtSPforRefFile.squaredDistanceSum(i) = gammaParam.squaredDistanceSum;
        infoAtSPforRefFile.HSAVERS_SPlat(i) = SP_lla(i,1);
        infoAtSPforRefFile.HSAVERS_SPlon(i) = SP_lla(i,2);
        infoAtSPforRefFile.HSAVERS_SPalt_5km(i) = meanSurfaceElevation;
        infoAtSPforRefFile.HSAVERS_SPalt_origDEM(i) = SPelevFromDEM(i);
        if indexTxFrequency==1
            infoAtSPforRefFile.RxGainLR1atSP(i) = gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency);
            infoAtSPforRefFile.RxGainRR1atSP(i) = gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency);
            infoAtSPforRefFile.Tx1GainAtSP(i)   = TRAp.TxGainAtSP(indexTxFrequency);
            infoAtSPforRefFile.Tx1GainAtRx(i)   = TxGainAtRx(i,indexTxFrequency);
            delayAtSPReferencedToDirect_1_microsec(i) = DDMparam.delayAtSPReferencedToDirect_microsec(indexTxFrequency);
            absoluteDelayAtSP_1_microsec(i)           = DDMparam.absoluteDelayAtSP_microsec(indexTxFrequency);

            infoAtSPforRefFile.Yaw_TX(i)   = TRAp.Yaw_TX;
            infoAtSPforRefFile.Yaw_TX_Error(i)   = TRAp.Yaw_TX_error;

        else
            infoAtSPforRefFile.RxGainLR5atSP(i) = gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency);
            infoAtSPforRefFile.RxGainRR5atSP(i) = gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency);
            infoAtSPforRefFile.Tx5GainAtSP(i) = TRAp.TxGainAtSP(indexTxFrequency);
            infoAtSPforRefFile.Tx5GainAtRx(i) = TxGainAtRx(i,indexTxFrequency);
            delayAtSPReferencedToDirect_5_microsec(i) = DDMparam.delayAtSPReferencedToDirect_microsec(indexTxFrequency);
            absoluteDelayAtSP_5_microsec(i) = DDMparam.absoluteDelayAtSP_microsec(indexTxFrequency);

            %infoAtSPforRefFile.Yaw_TX(i)   = TRAp.Yaw_TX;

        end
        %% ----------- 
   
        fclose('all');
    end
end

%% Coherent Channel Interpolation
%

lat_sp = SP_lla(:,1); % Lat of Specular Point
lon_sp = SP_lla(:,2); % Lon of Specular Point

M      = geoSYSp.incoIntTime./geoSYSp.cohIntTime;
    
%waveL  = TRAp.waveL;

for indexTxFrequency=1:length(TRAp.TxFrequency)
    %it evaluate the coherent component only in trackwise mode
    if strcmp(geoSYSp.Mode,'monte carlo')  == 1 && tag.AllSP==1
        if Nruns > 1
            
            %Original Version
            %AdHocFactor   = 1*10^8;%1250000/(4*pi);

            slowtime      = 1:size(AcohLR,1);
            fasttime      = (1:(M*size(AcohLR,1)))./M;

             AcohAbsFine   = interp1(slowtime, abs(AcohLR(:,indexTxFrequency)),  fasttime, 'spline');
            AincohAbsFine = interp1(slowtime, abs(AdiffLR(:,indexTxFrequency)), fasttime, 'spline');
            DistFine      = interp1(slowtime, rho_TRsp(:,indexTxFrequency), fasttime, 'spline');

            %New version - June 2023
            AdHocFactor   = 1*10^10./(AcohAbsFine.^2 + 1); 

            lat_spFine    = interp1(slowtime, lat_sp, fasttime, 'spline');
            lon_spFine    = interp1(slowtime, lon_sp, fasttime, 'spline');

        	% Coherent complex field at high-sampling rate:
        	AcohFine = AcohAbsFine.*exp((1j*2*pi/TRAp.waveL(indexTxFrequency)).*DistFine);

        	% Incoherent complex field at high-sampling rate,
        	% 1) generate an array of random values that follow an exponential
        	%    distribution, of same length as the fine time series

            power_exponential = zeros(1, length(AincohAbsFine));

            for k = 1:length(AincohAbsFine)
                power_exponential(k) = random('Exponential',(AincohAbsFine(k))^2);
            end

        	AincohFine = AdHocFactor .* sqrt(power_exponential) .* exp((1j*2*pi) .* rand(1, length(AincohAbsFine)));
        	Atotal     = AcohFine + AdHocFactor.*AincohFine;
            outputSimulations.lat_spFine = lat_spFine;
            outputSimulations.lat_spFine = lon_spFine;
            
            % %Coherent Channel Savings
            % labelsave = 'T1_500m_noF';
            % 
            % save(['AcohFine_',labelsave, num2str(i),'.mat'],'AcohFine')
            % save(['AincohFine_',labelsave, num2str(i),'.mat'],'AincohFine')
            % 
            % save(['lat_spFine_',labelsave, num2str(i),'.mat'],'lat_spFine')
            % save(['lon_spFine_',labelsave, num2str(i),'.mat'],'lon_spFine')
            % 
            % save(['lat_sp_',labelsave,'.mat'],'lat_sp')
            % save(['lon_sp_',labelsave,'.mat'],'lon_sp')     
            % 
            % save(['DistFin_',labelsave,'.mat'],'DistFine')   

            if indexTxFrequency==1
                outputSimulations.Atotal_1     = Atotal;
                outputSimulations.Atotal_5     = zeros(M*size(AcohLR,1),1);
            else
                outputSimulations.Atotal_5     = Atotal;
            end

         end
     else
         
         outputSimulations.Atotal_1     = zeros(M*size(AcohLR,1),1);
         outputSimulations.Atotal_5     = zeros(M*size(AcohLR,1),1);
         
         outputSimulations.lat_spFine = zeros(M*size(AcohLR,1),1);
         outputSimulations.lat_spFine = zeros(M*size(AcohLR,1),1);
     end
end

       %mean(abs(AincohFine(1:250)).^2)
       %std(abs(AincohFine(1:250)).^2)
 
       %figure, plot(AincohFine,'o'), xlabel('Re'), ylabel('Im')
       %figure, plot(AcohFine,'o'), xlabel('Re'), ylabel('Im')

       %figure, plot(10^11*abs(AincohFine)), xlabel('index'), ylabel('|inco|')
       %figure, plot(abs(AcohFine)), xlabel('index'), ylabel('|coh|')

       %figure, plot(abs(10^12.*AincohFine + AcohFine)), xlabel('index'), ylabel('|coh|')
       %figure, plot(angle(10^12.*AincohFine + AcohFine)), xlabel('index'), ylabel('|coh|')

       %AincohFine_c = AincohFine.*exp((-1j*2*pi/TRAp.waveL).*DistFine);
       %AcohFine_c   = AcohFine.*exp((-1j*2*pi/TRAp.waveL).*DistFine);

      %figure, plot(angle(10^0.*AincohFine_c + AcohFine_c)), xlabel('index'), ylabel('angle(coh)')

       %figure, plot(angle(Atotal)), xlabel('index'), ylabel('angle(Atot)')
       %figure, plot(abs(Atotal)), xlabel('index'), ylabel('|Atot|')

       %save('lat_spFine.mat', 'lat_spFine')
       %save('lon_spFine.mat', 'lon_spFine')
       %save('AincohFine.mat', 'AincohFine')
       %save('AcohFine.mat', 'AcohFine')
       %save('DistFine.mat', 'DistFine')
       %save('waveL.mat','waveL')    



%% SAVE THE OUTPUTS IN A STRUCTURE, L1/E1

outputSimulations.max_sigPower_noNoise_LR1 = max_sigPower_noNoise_NoError_LR(:,1);
outputSimulations.max_sigPower_noNoise_RR1 = max_sigPower_noNoise_NoError_RR(:,1);
outputSimulations.max_sigPower_Noise_LR1   = max_sigPower_Noise_LR(:,1);
outputSimulations.max_sigPower_Noise_RR1   = max_sigPower_Noise_RR(:,1);

outputSimulations.sigPower_noNoise_LR1_lin = sigPower_noNoise_NoError_LR_lin(:,:,:,1);
outputSimulations.sigPower_noNoise_RR1_lin = sigPower_noNoise_NoError_RR_lin(:,:,:,1);
outputSimulations.sigPower_Noise_LR1_lin   = sigPower_Noise_LR_lin(:,:,:,1);
outputSimulations.sigPower_Noise_RR1_lin   = sigPower_Noise_RR_lin(:,:,:,1);

if tag.SNR==1
    outputSimulations.SNR_LR1 = SNR_LR(:,1);
    outputSimulations.SNR_RR1 = SNR_RR(:,1);
%     outputSimulations.geobSNR_LR1 = geobSNR_LR(:,1);
%     outputSimulations.geobSNR_RR1 = geobSNR_RR(:,1);    
end

outputSimulations.delayAtSPReferencedToDirect_1_microsec = delayAtSPReferencedToDirect_1_microsec;
outputSimulations.absoluteDelayAtSP_1_microsec           = absoluteDelayAtSP_1_microsec;

%%%power black body nadir DDM, LR
outputSimulations.NLblackbody_DOWN_LR_1     = NLbb_DOWN_LR(:,:,:,1);
outputSimulations.NLblackbody_DOWNmean_LR_1 = squeeze(mean(mean(outputSimulations.NLblackbody_DOWN_LR_1)));
outputSimulations.NLblackbody_DOWNmax_LR_1  = squeeze(max(max(outputSimulations.NLblackbody_DOWN_LR_1)));

%%%power black body nadir DDM, RR
outputSimulations.NLblackbody_DOWN_RR_1     = NLbb_DOWN_RR(:,:,:,1);
outputSimulations.NLblackbody_DOWNmean_RR_1 = squeeze(mean(mean(outputSimulations.NLblackbody_DOWN_RR_1)));
outputSimulations.NLblackbody_DOWNmax_RR_1  = squeeze(max(max(outputSimulations.NLblackbody_DOWN_RR_1)));

%%%power black body zenith
outputSimulations.NLblackbody_UP_1     = NLbb_UP(:,:,:,1);
outputSimulations.NLblackbody_UPmean_1 = squeeze(mean(mean(outputSimulations.NLblackbody_UP_1)));
outputSimulations.NLblackbody_UPmax_1  = squeeze(max(max(outputSimulations.NLblackbody_UP_1)));

%%%power antenna UP
outputSimulations.LR1_UP_Power = LR_UP_Power(:,1);
outputSimulations.RR1_UP_Power = RR_UP_Power(:,1);

%%%power antenna UP Flu
outputSimulations.LR1_UP_Power_Fluc = LR_UP_Power_Flu(:,1);
outputSimulations.RR1_UP_Power_Fluc = RR_UP_Power_Flu(:,1);
outputSimulations.PN_star_1 = PN_star(:,1);

outputSimulations.max_reflectivity_noNoise_LR1_dB=max_reflectivity_noNoise_NoError_LR_tot_dB(:,1);
outputSimulations.max_reflectivity_noNoise_RR1_dB=max_reflectivity_noNoise_NoError_RR_tot_dB(:,1);
outputSimulations.max_reflectivity_Noisy_LR1_dB=max_reflectivity_Noisy_LR_tot_dB(:,1);
outputSimulations.max_reflectivity_Noisy_RR1_dB=max_reflectivity_Noisy_RR_tot_dB(:,1);
outputSimulations.max_geobReflec_noNoise_LR1_dB = max_geobReflec_noNoise_NoError_LR_dB(:,1);
outputSimulations.max_geobReflec_noNoise_RR1_dB = max_geobReflec_noNoise_NoError_RR_dB(:,1);
outputSimulations.max_geobPower_Noisy_LR1 = max_geobPower_Noisy_LR(:,1);
outputSimulations.max_geobPower_Noisy_RR1 = max_geobPower_Noisy_RR(:,1);

%% SAVE THE OUTPUTS IN A STRUCTURE, L5/E5
outputSimulations.max_sigPower_noNoise_LR5 = max_sigPower_noNoise_NoError_LR(:,2);
outputSimulations.max_sigPower_noNoise_RR5 = max_sigPower_noNoise_NoError_RR(:,2);
outputSimulations.max_sigPower_Noise_LR5   = max_sigPower_Noise_LR(:,2);
outputSimulations.max_sigPower_Noise_RR5   = max_sigPower_Noise_RR(:,2);

outputSimulations.sigPower_noNoise_LR5_lin = sigPower_noNoise_NoError_LR_lin(:,:,:,2);
outputSimulations.sigPower_noNoise_RR5_lin = sigPower_noNoise_NoError_RR_lin(:,:,:,2);
outputSimulations.sigPower_Noise_LR5_lin   = sigPower_Noise_LR_lin(:,:,:,2);
outputSimulations.sigPower_Noise_RR5_lin   = sigPower_Noise_RR_lin(:,:,:,2);

if tag.SNR==1
    outputSimulations.SNR_LR5 = SNR_LR(:,2);
    outputSimulations.SNR_RR5 = SNR_RR(:,2);
%     outputSimulations.geobSNR_LR5 = geobSNR_LR(:,2);
%     outputSimulations.geobSNR_RR5 = geobSNR_RR(:,2);    
end

outputSimulations.delayAtSPReferencedToDirect_5_microsec = delayAtSPReferencedToDirect_5_microsec;
outputSimulations.absoluteDelayAtSP_5_microsec = absoluteDelayAtSP_5_microsec;

%%%power black body nadir DDM, LR
outputSimulations.NLblackbody_DOWN_LR_5     = NLbb_DOWN_LR(:,:,:,2);
outputSimulations.NLblackbody_DOWNmean_LR_5 = squeeze(mean(mean(outputSimulations.NLblackbody_DOWN_LR_5)));
outputSimulations.NLblackbody_DOWNmax_LR_5  = squeeze(max(max(outputSimulations.NLblackbody_DOWN_LR_5)));

%%%power black body nadir DDM, RR
outputSimulations.NLblackbody_DOWN_RR_5     = NLbb_DOWN_RR(:,:,:,2);
outputSimulations.NLblackbody_DOWNmean_RR_5 = squeeze(mean(mean(outputSimulations.NLblackbody_DOWN_RR_5)));
outputSimulations.NLblackbody_DOWNmax_RR_5  = squeeze(max(max(outputSimulations.NLblackbody_DOWN_RR_5)));

%%%power black body zenith
outputSimulations.NLblackbody_UP_5     = NLbb_UP(:,:,:,2);
outputSimulations.NLblackbody_UPmean_5 = squeeze(mean(mean(outputSimulations.NLblackbody_UP_5)));
outputSimulations.NLblackbody_UPmax_5  = squeeze(max(max(outputSimulations.NLblackbody_UP_5)));

%%%power antenna UP
outputSimulations.LR5_UP_Power = LR_UP_Power(:,2);
outputSimulations.RR5_UP_Power = RR_UP_Power(:,2);

%%%power antenna UP Flu
outputSimulations.LR5_UP_Power_Fluc = LR_UP_Power_Flu(:,2);
outputSimulations.RR5_UP_Power_Fluc = RR_UP_Power_Flu(:,2);
outputSimulations.PN_star_5 = PN_star(:,2);

outputSimulations.max_reflectivity_noNoise_LR5_dB=max_reflectivity_noNoise_NoError_LR_tot_dB(:,2);
outputSimulations.max_reflectivity_noNoise_RR5_dB=max_reflectivity_noNoise_NoError_RR_tot_dB(:,2);
outputSimulations.max_reflectivity_Noisy_LR5_dB=max_reflectivity_Noisy_LR_tot_dB(:,2);
outputSimulations.max_reflectivity_Noisy_RR5_dB=max_reflectivity_Noisy_RR_tot_dB(:,2);
outputSimulations.max_geobReflec_noNoise_LR5_dB = max_geobReflec_noNoise_NoError_LR_dB(:,2);
outputSimulations.max_geobReflec_noNoise_RR5_dB = max_geobReflec_noNoise_NoError_RR_dB(:,2);
outputSimulations.max_geobPower_Noisy_LR5 = max_geobPower_Noisy_LR(:,2);
outputSimulations.max_geobPower_Noisy_RR5 = max_geobPower_Noisy_RR(:,2);


%% =============== END of Electromagnetic Simulation =====================

    end