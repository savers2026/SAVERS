
%% ************************************************************************
% MODULE NAME:      selectAndReadGNSSRdata_newOrbitFile.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, N. Pierdicca, A. Musunuri
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      Sep 2019 by Laura Dente
% Main updates:          2020 – 2022 by L. Dente, N. Pierdicca, Aneesha Musunuri
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module reads the orbit file and extracts the information related to
% the specular points selected by the user, i.e. a specific SP, or a track,
% or all (or a selection of) SPs in a specific time range.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [Tx_selected, VTx_selected, Rx_selected, VRx_selected,filenameFigDEMandDDM,...
            SPlatlonAlt_selectedSP,incAngle_selectedSP,PRN_selected,...
            GNSS_Week_selected,GNSS_SoW_selected,trackID_selected,SVN_selected,goToSimulation,OrbitWithReflectivity,ReflectivitydB_selected,tag,Rx_roll_selected,Rx_pitch_selected,Rx_yaw_selected,Rx_Attitude_Variation_RPY]= ...
            selectAndReadGNSSRdata_newOrbitFile(filenameData,DEMinfo,filenameCoverMap,outputDirectory, mode1, geoSYSp,Sampling,Sampling_all_Text,tag)

switch mode1
    case 'Single Point'
        flagMode=1;
    case 'Trackwise'
        flagMode=2;
    case 'Global'
        flagMode=3;
    case 'Monte Carlo'
        flagMode=4;
end

%to identify orbit file with Reflectivity
infoOrbitFile=ncinfo(filenameData);
allVarNames = {infoOrbitFile.Variables.Name};
nameMatch = strncmpi(allVarNames,'ReflectivitydB',14);
OrbitWithReflectivity=sum(nameMatch);

%read the parameters repoted in the GNSSR netcdf file

%case ALL DAYS
if strcmp(Sampling_all_Text{1},'All Days')

    TxID    = ncread(filenameData,'GNSS_System'); %'E' is Galileo 'G' is GPS
    PRN     = ncread(filenameData,'PRN');
    SVN     = ncread(filenameData,'SVN');
    trackID = ncread(filenameData,'Track_ID');
    SPlat   = ncread(filenameData,'Lat_SP');
    SPlon   = ncread(filenameData,'Lon_SP');
    SPAlt    = ncread(filenameData,'Alt_SP');
    incAngle =ncread(filenameData,'incidence');
    numChannels=ncread(filenameData,'N_channels'); %4 or 10 channels
    %time_temp=ncread(filenameData,'OBS_Time');
    Tx_x=ncread(filenameData,'x_tx');
    Tx_y=ncread(filenameData,'y_tx');
    Tx_z=ncread(filenameData,'z_tx');
    Tx_Vx=ncread(filenameData,'vx_tx');
    Tx_Vy=ncread(filenameData,'vy_tx');
    Tx_Vz=ncread(filenameData,'vz_tx');
    Rx_x_temp=ncread(filenameData,'x_rcv');
    Rx_y_temp=ncread(filenameData,'y_rcv');
    Rx_z_temp=ncread(filenameData,'z_rcv');
    Rx_Vx_temp=ncread(filenameData,'vx_rcv');
    Rx_Vy_temp=ncread(filenameData,'vy_rcv');
    Rx_Vz_temp=ncread(filenameData,'vz_rcv');
    Rx_lat_temp = ncread(filenameData,'Lat_rcv'); %Amirreza
    Rx_lon_temp = ncread(filenameData,'Lon_rcv'); %Amirreza

    
    % Variation Rx attitude by Amirreza
    Rx_Attitude_Variation_RPY = [0 0 0];  % [roll pitch yaw]
    if ismember('roll_rcv', allVarNames)
            Rx_Attitude_Variation_RPY(1) = 1;
    Rx_roll_temp = ncread(filenameData, 'roll_rcv');
    else
        Rx_roll_temp = zeros(size(Rx_lat_temp));
        Rx_Attitude_Variation_RPY(1) = 0;
    end
    
    if ismember('pitch_rcv', allVarNames)
        Rx_pitch_temp = ncread(filenameData, 'pitch_rcv');
            Rx_Attitude_Variation_RPY(2) = 1;
    else
        Rx_pitch_temp = zeros(size(Rx_lat_temp));
            Rx_Attitude_Variation_RPY(2) = 0;
    end
    
    if ismember('yaw_rcv', allVarNames)
        Rx_yaw_temp = ncread(filenameData, 'yaw_rcv');
            Rx_Attitude_Variation_RPY(3) = 1;
    else
        Rx_yaw_temp = zeros(size(Rx_lat_temp));
            Rx_Attitude_Variation_RPY(3) = 0;
    end
    % Variation Rx attitude by Amirreza

    GNSS_Week_temp=ncread(filenameData,'GPS_Week');
    GNSS_SoW_temp=ncread(filenameData,'GPS_SoW');
%     GnssBlock=ncread(filenameData,'GnssBlock');
    if OrbitWithReflectivity
        ReflectivitydB=ncread(filenameData, 'ReflectivitydB'); 
    end

%case TIME RANGE in input 
else
    Init_day = Sampling{1};
    Last_day = Sampling{2};

    GNSS_Week_temp=ncread(filenameData,'GPS_Week');
    GNSS_SoW_temp=ncread(filenameData,'GPS_SoW');
%     GnssBlock=ncread(filenameData,'GnssBlock');
    Orbit_days=fix(GNSS_Week_temp(:)*7+GNSS_SoW_temp(:)/86400-GNSS_Week_temp(1)*7-GNSS_SoW_temp(1)/86400)+1;
    if Init_day<Orbit_days(1) || Orbit_days(end)<Init_day
        Init_day=1;
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'WARNING: The selected initial day exceeds the orbit length.';...
                            ['It will be set to ' num2str(Init_day)]},CreateStruct))
    end

    if Last_day<Orbit_days(1) || Orbit_days(end)<Last_day
        Last_day=Orbit_days(end);
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'WARNING: The selected last day exceeds the orbit length.';...
                            ['It will be set to ' num2str(Last_day)]},CreateStruct))
    end

    Init_day=find(Orbit_days>=Init_day, 1) ;
    Last_day=find(Orbit_days<=Last_day, 1, 'last') ;
    Num_of_elements=Last_day-Init_day+1 ;
    GNSS_Week_temp=GNSS_Week_temp(Init_day:Last_day);
    GNSS_SoW_temp=GNSS_SoW_temp(Init_day:Last_day);
    % read the orbit file variables in the selected days
    numChannels=ncread(filenameData,'N_channels'); %4 or 10 channels
    TxID=ncread(filenameData,'GNSS_System',[Init_day 1], [Num_of_elements numChannels]); %'E' is Galileo 'G' is GPS
    PRN=ncread(filenameData,'PRN',[Init_day 1], [Num_of_elements numChannels]);
    SVN=ncread(filenameData,'SVN',[Init_day 1], [Num_of_elements numChannels]);
%     GnssBlock=ncread(filenameData,'GnssBlock',[Init_day 1], [Num_of_elements numChannels]);
    trackID=ncread(filenameData,'Track_ID',[Init_day 1], [Num_of_elements numChannels]);
    SPlat=ncread(filenameData,'Lat_SP',[Init_day 1], [Num_of_elements numChannels]);
    SPlon=ncread(filenameData,'Lon_SP',[Init_day 1], [Num_of_elements numChannels]);
    SPAlt=ncread(filenameData,'Alt_SP',[Init_day 1], [Num_of_elements numChannels]);
    incAngle=ncread(filenameData,'incidence',[Init_day 1], [Num_of_elements numChannels]);
    %time_temp=ncread(filenameData,'OBS_Time');
    Tx_x=ncread(filenameData,'x_tx',[Init_day 1], [Num_of_elements numChannels]);
    Tx_y=ncread(filenameData,'y_tx',[Init_day 1], [Num_of_elements numChannels]);
    Tx_z=ncread(filenameData,'z_tx',[Init_day 1], [Num_of_elements numChannels]);
    Tx_Vx=ncread(filenameData,'vx_tx',[Init_day 1], [Num_of_elements numChannels]);
    Tx_Vy=ncread(filenameData,'vy_tx',[Init_day 1], [Num_of_elements numChannels]);
    Tx_Vz=ncread(filenameData,'vz_tx',[Init_day 1], [Num_of_elements numChannels]);
    Rx_x_temp=ncread(filenameData,'x_rcv',Init_day, Num_of_elements);
    Rx_y_temp=ncread(filenameData,'y_rcv',Init_day, Num_of_elements);
    Rx_z_temp=ncread(filenameData,'z_rcv',Init_day, Num_of_elements);
    Rx_Vx_temp=ncread(filenameData,'vx_rcv',Init_day, Num_of_elements);
    Rx_Vy_temp=ncread(filenameData,'vy_rcv',Init_day, Num_of_elements);
    Rx_Vz_temp=ncread(filenameData,'vz_rcv',Init_day, Num_of_elements);
%     Rx_lat_temp=ncread(filenameData,'Lat_rcv',Init_day, Num_of_elements); %Amirreza
%     Rx_lon_temp=ncread(filenameData,'Lon_rcv',Init_day, Num_of_elements); %Amirreza
    Rx_lat_temp=ncread(filenameData,'Lat_rcv',Init_day, Num_of_elements); %Amirreza
    Rx_lon_temp=ncread(filenameData,'Lon_rcv',Init_day, Num_of_elements); %Amirreza

    % Variation Rx attitude by Amirreza

    if ismember('roll_rcv', allVarNames)
    Rx_roll_temp = ncread(filenameData, 'roll_rcv', Init_day, Num_of_elements);
        Rx_Attitude_Variation_RPY(1) = 1;
    else
        Rx_roll_temp = zeros(size(Rx_lat_temp));
            Rx_Attitude_Variation_RPY(1) = 0;
    end
    
    if ismember('pitch_rcv', allVarNames)
        Rx_pitch_temp = ncread(filenameData, 'pitch_rcv', Init_day, Num_of_elements);
            Rx_Attitude_Variation_RPY(2) = 1;
    else
        Rx_pitch_temp = zeros(size(Rx_lat_temp));
            Rx_Attitude_Variation_RPY(2) = 0;
    end
    
    if ismember('yaw_rcv', allVarNames)
        Rx_yaw_temp = ncread(filenameData, 'yaw_rcv', Init_day, Num_of_elements);
            Rx_Attitude_Variation_RPY(3) = 1;
    else
        Rx_yaw_temp = zeros(size(Rx_lat_temp));
            Rx_Attitude_Variation_RPY(3) = 0;
    end
        
    % Variation Rx attitude by Amirreza

    if OrbitWithReflectivity
        ReflectivitydB=ncread(filenameData, 'ReflectivitydB',[Init_day 1], [Num_of_elements numChannels]);
    end
end

GNSS_Week=repmat(GNSS_Week_temp,1,numChannels)';
GNSS_SoW=repmat(GNSS_SoW_temp,1,numChannels)';
Rx_x=repmat(Rx_x_temp,1,numChannels)';
Rx_y=repmat(Rx_y_temp,1,numChannels)';
Rx_z=repmat(Rx_z_temp,1,numChannels)';
Rx_Vx=repmat(Rx_Vx_temp,1,numChannels)';
Rx_Vy=repmat(Rx_Vy_temp,1,numChannels)';
Rx_Vz=repmat(Rx_Vz_temp,1,numChannels)';
Rx_lat = repmat(Rx_lat_temp,1,numChannels)';
Rx_lon = repmat(Rx_lon_temp,1,numChannels)';

Rx_roll = repmat(Rx_roll_temp,1,numChannels)'; %Amirreza
Rx_pitch = repmat(Rx_pitch_temp,1,numChannels)'; %Amirreza
Rx_yaw = repmat(Rx_yaw_temp,1,numChannels)'; %Amirreza


%resize the variables for the area of interest and for the selected
%trasmitter (GPS or Galileo)
if geoSYSp.indexTxSignal==0 || geoSYSp.indexTxSignal==1 %GPS
    indexSPLatLonSelectedSiteAndTx= find(SPlat>=DEMinfo.SPcoodinatesLimits(1) & SPlat<=DEMinfo.SPcoodinatesLimits(2) &...
            SPlon>=DEMinfo.SPcoodinatesLimits(3) & SPlon<=DEMinfo.SPcoodinatesLimits(4) & TxID=='G');
else %Galileo
    indexSPLatLonSelectedSiteAndTx= find(SPlat>=DEMinfo.SPcoodinatesLimits(1) & SPlat<=DEMinfo.SPcoodinatesLimits(2) &...
            SPlon>=DEMinfo.SPcoodinatesLimits(3) & SPlon<=DEMinfo.SPcoodinatesLimits(4) & TxID=='E');
end

SPlat_selectedSiteAndTx=SPlat(indexSPLatLonSelectedSiteAndTx);

%% Amirreza
% if size(TxID) ~= size(Rx_lon)
%     TxID = TxID';
% end
% 
% if geoSYSp.indexTxSignal==0 || geoSYSp.indexTxSignal==1 %GPS
%     indexRxLatLonSelectedSiteAndTx= find(Rx_lat>=DEMinfo.SPcoodinatesLimits(1) & Rx_lat<=DEMinfo.SPcoodinatesLimits(2) &...
%             Rx_lon>=DEMinfo.SPcoodinatesLimits(3) & Rx_lon<=DEMinfo.SPcoodinatesLimits(4) & TxID=='G');
% else %Galileo
%     indexRxLatLonSelectedSiteAndTx= find(Rx_lat>=DEMinfo.SPcoodinatesLimits(1) & Rx_lat<=DEMinfo.SPcoodinatesLimits(2) &...
%             Rx_lon>=DEMinfo.SPcoodinatesLimits(3) & Rx_lon<=DEMinfo.SPcoodinatesLimits(4) & TxID=='E');
% end
% 
% Rxlat_selectedSiteAndTx=Rx_lat(indexRxLatLonSelectedSiteAndTx);

%% Amirreza

if isempty(SPlat_selectedSiteAndTx)
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    uiwait(msgbox({'Input error: The orbit file contains zero specular points in the selected time slot, over the selected site and for the selected transmitter.';...
                'The simulation will be stopped.';...
                'Please run again the simulation selecting either a different time slot or a different site or a different transmitter'},'INPUT ERROR MESSAGE',CreateStruct))
    Tx_selected=0; VTx_selected=0;Rx_selected=0;VRx_selected=0;SPlatlonAlt_selectedSP=[0,0,0];incAngle_selectedSP=0;
    PRN_selected=0;GNSS_Week_selected=0;GNSS_SoW_selected=0;trackID_selected=0;SVN_selected=0;goToSimulation=0;
    OrbitWithReflectivity=0;ReflectivitydB_selected=0;Rx_roll_selected=0;Rx_pitch_selected=0;Rx_yaw_selected=0;
    filenameFigDEMandDDM='';
    close all
    return
end


Rxlat_selectedSiteAndTx=Rx_lat(indexSPLatLonSelectedSiteAndTx);
Rxlon_selectedSiteAndTx=Rx_lon(indexSPLatLonSelectedSiteAndTx);
SPlon_selectedSiteAndTx=SPlon(indexSPLatLonSelectedSiteAndTx);
SPAlt_selectedSiteAndTx=SPAlt(indexSPLatLonSelectedSiteAndTx);
incAngle_selectedSiteAndTx=incAngle(indexSPLatLonSelectedSiteAndTx);
PRN_selectedSiteAndTx=PRN(indexSPLatLonSelectedSiteAndTx);
SVN_selectedSiteAndTx=SVN(indexSPLatLonSelectedSiteAndTx);
trackID_selectedSiteAndTx=trackID(indexSPLatLonSelectedSiteAndTx);
GNSS_Week_selectedSiteAndTx=GNSS_Week(indexSPLatLonSelectedSiteAndTx);
GNSS_SoW_selectedSiteAndTx=GNSS_SoW(indexSPLatLonSelectedSiteAndTx);
% GnssBlock_selectedSiteAndTx=GnssBlock(indexSPLatLonSelectedSiteAndTx);
TxID_selectedSiteAndTx=TxID(indexSPLatLonSelectedSiteAndTx);
Tx_x_selectedSiteAndTx=Tx_x(indexSPLatLonSelectedSiteAndTx);
Tx_y_selectedSiteAndTx=Tx_y(indexSPLatLonSelectedSiteAndTx);
Tx_z_selectedSiteAndTx=Tx_z(indexSPLatLonSelectedSiteAndTx);
Tx_Vx_selectedSiteAndTx=Tx_Vx(indexSPLatLonSelectedSiteAndTx);
Tx_Vy_selectedSiteAndTx=Tx_Vy(indexSPLatLonSelectedSiteAndTx);
Tx_Vz_selectedSiteAndTx=Tx_Vz(indexSPLatLonSelectedSiteAndTx);
Rx_x_selectedSiteAndTx=Rx_x(indexSPLatLonSelectedSiteAndTx);
Rx_y_selectedSiteAndTx=Rx_y(indexSPLatLonSelectedSiteAndTx);
Rx_z_selectedSiteAndTx=Rx_z(indexSPLatLonSelectedSiteAndTx);
Rx_Vx_selectedSiteAndTx=Rx_Vx(indexSPLatLonSelectedSiteAndTx);
Rx_Vy_selectedSiteAndTx=Rx_Vy(indexSPLatLonSelectedSiteAndTx);
Rx_Vz_selectedSiteAndTx=Rx_Vz(indexSPLatLonSelectedSiteAndTx);

Rx_roll_selectedSiteAndTx=Rx_roll(indexSPLatLonSelectedSiteAndTx); %Amirreza
Rx_pitch_selectedSiteAndTx=Rx_pitch(indexSPLatLonSelectedSiteAndTx); %Amirreza
Rx_yaw_selectedSiteAndTx=Rx_yaw(indexSPLatLonSelectedSiteAndTx); %Amirreza

if OrbitWithReflectivity 
    ReflectivitydB_selectedSiteAndTx=ReflectivitydB(indexSPLatLonSelectedSiteAndTx);
end
             
%% ------SELECT THE START POINT OF THE SP TRACK TO BE SIMULATED------------

if flagMode==1 || flagMode==2 % Single point or Trackwise

    if tag.AllSP==0
        tag.AllSP=1;
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Warning: The subsampling of the SPs number is not possible in Single Point and Trackwise mode.';...
                            'All SPs available in the selected time range will be shown.'},CreateStruct))
    end

    if DEMinfo.flatsurface==0
        %open and read the DEM of the area of interest       
        [DEM_lla,~,~] = readDEMforPlot(DEMinfo.filenameDEM,DEMinfo.SPcoodinatesLimits);
        %show DEM in lat-lon and overlap the SP tracks
        figure('Name', 'DEM', 'NumberTitle','off','OuterPosition', [150 550 600 500],'Units','normalized'); %[50 300 600 500]
        mapshow(DEM_lla(:,:,2),DEM_lla(:,:,1),DEM_lla(:,:,3), 'DisplayType', 'texturemap')
        axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
        cbarDEM=colorbar;
        cbarDEM.Label.String = 'Surface elevation (m)';
        title('DEM - SP white circle') 
        xlabel('Longitude (deg)')
        ylabel('Latitude (deg)')
        hold on
        scatter(SPlon_selectedSiteAndTx,SPlat_selectedSiteAndTx,25,[1,1,1])
    end
    
    %% Amirreza

    % Showing RX path (Track) and SP around it
    figure('Name', 'Rx Track', 'NumberTitle', 'off', 'OuterPosition', [750 550 300 500], 'Units', 'normalized');
    gx = geoaxes;
    geobasemap satellite;
    gx.TickLabelFormat = '-dd';
    zoom on;
    hold on;
    
    unique_track_ids = unique(trackID_selectedSiteAndTx);
    num_tracks = length(unique_track_ids);
    cmap = lines(num_tracks);
    
    darker_factor = 0.2;  % Reduce brightness for Rx track (0 = black)
    lighter_factor = 0.7; % Increase brightness for SP points (1 = white)
    
    color_map_sp = cmap * (1 - lighter_factor) + lighter_factor;  
    color_map_rx = cmap * (1 - darker_factor);
    
    % Plot SP track
    for i = 1:num_tracks
        idx = trackID_selectedSiteAndTx == unique_track_ids(i);
        geoscatter(SPlat_selectedSiteAndTx(idx), SPlon_selectedSiteAndTx(idx), 30, ...
            'MarkerEdgeColor', color_map_sp(i,:), 'MarkerFaceColor', color_map_sp(i,:), ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7, 'Parent', gx);
    end
    
    % Plot Rx track
    for i = 1:num_tracks
        idx = trackID_selectedSiteAndTx == unique_track_ids(i);
        geoscatter(Rxlat_selectedSiteAndTx(idx), Rxlon_selectedSiteAndTx(idx), 60, ...
            'MarkerEdgeColor', color_map_rx(i,:), 'MarkerFaceColor', color_map_rx(i,:), ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1, 'Parent', gx);
    end
    
    title('Rx Tracks (Darker) Vs SP (Lighter)');


    %% Amirreza

    %plot SP tracks
    %the colour of each SP depends on the incidence angle
    figure('Name', 'Incidence angle', 'NumberTitle','off','OuterPosition', [1050 550 300 500],'Units','normalized');%[900 50 600 500]
    scatter(SPlon_selectedSiteAndTx,SPlat_selectedSiteAndTx,25,incAngle_selectedSiteAndTx,'filled')
    axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
    cbarIncAng=colorbar;
    cbarIncAng.Label.String = 'Incidence angle (°)';
    title('SP tracks vs Incidence angle')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')

    %open and read the land cover map of the area of interest
    SPlatlon_selectedSite(:,1)=SPlat_selectedSiteAndTx;
    SPlatlon_selectedSite(:,2)=SPlon_selectedSiteAndTx;
    [coverMap,~]=readCoverMapCCI(filenameCoverMap,SPlatlon_selectedSite,1.0);
    %plot the SP tracks over the land cover map of the area of interest
    mapOfColorCover=[0 0 1;0 1 0;1 0 1;0 1 1;1 0 0;1 1 0;0 0 0];
    figure('Name', 'Land cover map', 'NumberTitle','off','OuterPosition', [150 50 600 500],'Units','normalized');%[50 50 700 510]
    subplot('Position', [0.1 0.12 0.65 0.82]);
    mapshow(coverMap(:,:,2),coverMap(:,:,1),coverMap(:,:,3), 'DisplayType', 'texturemap')
    axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
    colormap(mapOfColorCover)
    colorbar('Ticks',[0,1,2,3,4,5,6],...
                'TickLabels',{'0 water bodies','1 broadleaved forest','2 niddleleaved forest','3 cropland/grassland/shrubland',...
                '4 bare areas','5 flooded forest','6 flooded shrubland'})
    caxis([0 6])
    title('Land cover map - SP white circle')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    hold on
    scatter(SPlon_selectedSiteAndTx,SPlat_selectedSiteAndTx,25,[1,1,1])        

    %plot SP tracks
    %the colour of each SP depends on the PRN
    %define colour map for PRN plot
    %{
    mapPRN=[1 1 0;0.5 0.5 0; 1 0 1;0.5 0 0.5;0 0.5 0.5; 0 1 1;0.5 1 0.5; 1 0 0;1 0.5 0.5; 0.5 0.5 1;...
            0 1 0;0.5 1 1; 1 0.5 1; 1 1 0.5; 0 0 1; 0.2 0.5 1; 0.5 1 0.2; 1 0.5 0.2; 0 0 0; 0.5 0.2 1;...
            0.7 0.5 1; 0.7 1 0.5; 1 0.7 0.5; 0.5 0.7 1; 0.5 1 0.7; 1 0.5 0.7; 0.2 0.5 0.7; 0.5 0.2 0.7; 0.5 0.7 0.2;...
            0.7 0.5 0.2; 0.7 0.2 0.5; 0.2 0.7 0.5; 0.3 0.9 0.1; 0.1 0.3 0.9; 0.9 0.1 0.3; 0.1 0.8 0.3];
    %}
    %% mapPRN color changed by Amirreza for GREAT PRNs and for better visibility
    mapPRN = lines(numel(unique(PRN_selectedSiteAndTx)));
    [~, ~, colorIdx] = unique(PRN_selectedSiteAndTx);
    
    figSPtracksPRN=figure('Name', 'PRN', 'NumberTitle','off','OuterPosition', [750 50 600 500],'Units','normalized');%[900 300 600 500]
    scatter(SPlon_selectedSiteAndTx,SPlat_selectedSiteAndTx,25,colorIdx,'filled')
    axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) ...
          DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
    colormap(mapPRN)
    caxis([0.5 numel(unique(PRN_selectedSiteAndTx))+0.5])
    cbarPRN=colorbar;
    cbarPRN.Ticks=1:numel(unique(PRN_selectedSiteAndTx));
    cbarPRN.TickLabels=string(unique(PRN_selectedSiteAndTx));
    cbarPRN.Label.String='PRN';
    title('SP coord vs PRN')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')


    % Note that the icon "data cursor" in the Figure bar menu should be
    % selected first
    datacursormode on
    dcm_obj = datacursormode(figSPtracksPRN);
    startMessageWindow=msgbox('Select starting specular point on PRN figure and then press OK (select first Data Cursor from the figure menu) ');
    uiwait(startMessageWindow);
    infoCursorStart_struct = getCursorInfo(dcm_obj);
    if isempty(infoCursorStart_struct)
        errordlg('No valid specular point selected','Input error');
        Tx_selected=0; VTx_selected=0;Rx_selected=0;VRx_selected=0;SPlatlonAlt_selectedSP=[0,0,0];incAngle_selectedSP=0;
        PRN_selected=0;GNSS_Week_selected=0;GNSS_SoW_selected=0;trackID_selected=0;SVN_selected=0;goToSimulation=0;
        OrbitWithReflectivity=0;ReflectivitydB_selected=0;Rx_roll_selected=0;Rx_pitch_selected=0;Rx_yaw_selected=0;
        filenameFigDEMandDDM='';
        close all
        return
    end
    %get latitude and data index of selected point
    SPlat_selectedSiteAndTx_all = SPlat_selectedSiteAndTx; %amir
    SPlon_selectedSiteAndTx_all = SPlon_selectedSiteAndTx; %amir

    latSP_start=infoCursorStart_struct.Position(1,2);
    lonSP_start=infoCursorStart_struct.Position(1,1); % checking SPlon Added by Amir for Grass/Great Project   
    indexSP_start=find(SPlat_selectedSiteAndTx_all==latSP_start &SPlon_selectedSiteAndTx_all==lonSP_start); % checking SPlon Added by Amir for Grass/Great Project   
    trackID_start=trackID_selectedSiteAndTx(indexSP_start);
    SPlat_selectedSiteAndTx_temp = SPlat_selectedSiteAndTx_all(trackID_selectedSiteAndTx == trackID_start); %amir
    SPlon_selectedSiteAndTx_temp = SPlon_selectedSiteAndTx_all(trackID_selectedSiteAndTx == trackID_start); %amir
    indexSP_start=find(SPlat_selectedSiteAndTx_temp==latSP_start &SPlon_selectedSiteAndTx_temp==lonSP_start); % checking SPlon Added by Amir for Grass/Great Project   

    %lonSP_start=infoCursorStart_struct.Position(1,1);
    %PRN_selectedSP=PRN_selectedSiteAndTx(indexSP_start);
    %SVN_selectedSP=SVN_selectedSiteAndTx(indexSP_start);
    %TxID_selectedSP=TxID_selectedSiteAndTx(indexSP_start); 
    %GNSS_Week_start=GNSS_Week_selectedSiteAndTx(indexSP_start);
    %GNSS_SoW_start=GNSS_SoW_selectedSiteAndTx(indexSP_start);
    %GnssBlock_selectedSP=GnssBlock_selectedSiteAndTx(indexSP_start);
    if length(indexSP_start)>1
        indexSP_start=indexSP_start(1);
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Warning: Multiple SPs with the same coordinates have been found.';...
                            'The first SP has been selected.'},CreateStruct))
    end


    %% select the end point of SP timeseries of interest
    switch mode1
        case 'Single Point'
            latSP_end=latSP_start;
            indexSP_end=indexSP_start;
            trackID_end = trackID_start;
        case 'Trackwise'
            trackCheckFlag=0;
            while trackCheckFlag<3        
                %select the end specular point
                datacursormode on
                dcm_obj = datacursormode(figSPtracksPRN);
                endMessageWindow=msgbox('Select ending specular point on the track previously selected in the PRN figure and then press OK (set Data Cursor from the figure menu) ');
                uiwait(endMessageWindow);
                infoCursorEnd_struct = getCursorInfo(dcm_obj);
                if isempty(infoCursorEnd_struct)
                    errordlg('No valid specular point selected','Input error');
                    Tx_selected=0; VTx_selected=0;Rx_selected=0;VRx_selected=0;SPlatlonAlt_selectedSP=[0,0,0];incAngle_selectedSP=0;
                    PRN_selected=0;GNSS_Week_selected=0;GNSS_SoW_selected=0;trackID_selected=0;SVN_selected=0;goToSimulation=0;
                    OrbitWithReflectivity=0;ReflectivitydB_selected=0;Rx_roll_selected=0;Rx_pitch_selected=0;Rx_yaw_selected=0;
                    filenameFigDEMandDDM='';
                    return
                end
                latSP_end=infoCursorEnd_struct.Position(1,2);   
                lonSP_end=infoCursorEnd_struct.Position(1,1);  % checking SPlon Added by Amir for Grass/Great Project   
                indexSP_end=find(SPlat_selectedSiteAndTx_all==latSP_end & SPlon_selectedSiteAndTx_all==lonSP_end); % checking SPlon Added by Amir for Grass/Great Project             
                trackID_end=trackID_selectedSiteAndTx(indexSP_end);

                SPlat_selectedSiteAndTx_temp = SPlat_selectedSiteAndTx_all(trackID_selectedSiteAndTx == trackID_end); %amir
                SPlon_selectedSiteAndTx_temp = SPlon_selectedSiteAndTx_all(trackID_selectedSiteAndTx == trackID_end); %amir
                indexSP_end=find(SPlat_selectedSiteAndTx_temp==latSP_end &SPlon_selectedSiteAndTx_temp==lonSP_end); % checking SPlon Added by Amir for Grass/Great Project   

                %lonSP_end=infoCursorEnd_struct.Position(1,1);
                %PRN_selectedSP_end=PRN_selectedSite(indexSP_end);
                if length(indexSP_end)>1
                    indexSP_end=indexSP_end(1);
                    CreateStruct.Interpreter = 'tex';
                    CreateStruct.WindowStyle = 'modal';
                    uiwait(msgbox({'Warning: Multiple SPs with the same coordinates have been found.';...
                            'The first SP has been selected.'},CreateStruct))
                end

                if trackID_end~=trackID_start
                    trackErrorMessageWindow=errordlg('The ending specular point is located on a different track than the starting specular point','Input error');
                    uiwait(trackErrorMessageWindow);
                    trackCheckFlag=trackCheckFlag+1;
                    continue
                else
                    trackCheckFlag=4;
                end
            end
            if trackCheckFlag~=4
                trackErrorMessageWindow2=errordlg('Wrong ending point selected. Please restart the simulation.','Input error');
                uiwait(trackErrorMessageWindow2);
                Tx_selected=0; VTx_selected=0;Rx_selected=0;VRx_selected=0;SPlatlonAlt_selectedSP=[0,0,0];incAngle_selectedSP=0;
                PRN_selected=0;GNSS_Week_selected=0;GNSS_SoW_selected=0;trackID_selected=0;SVN_selected=0;goToSimulation=0;
                OrbitWithReflectivity=0;ReflectivitydB_selected=0;Rx_roll_selected=0;Rx_pitch_selected=0;Rx_yaw_selected=0;
                filenameFigDEMandDDM='';
                close all
                return
            end
    end

    %% ----- extract all the info needed for the simulation for all selected SP -----
    %extract lat,lon,incidence angle of the selected series of specular points
    %Rx and Tx position and velocity

    %number selected SPs

         %amir
        trackID_selectedSiteAndTx_filtering = trackID_selectedSiteAndTx;
        % Rxlat_selectedSiteAndTx=Rx_lat(trackID_selectedSiteAndTx_filtering == trackID_end);
        % Rxlon_selectedSiteAndTx=Rx_lon(trackID_selectedSiteAndTx_filtering == trackID_end);
        SPlat_selectedSiteAndTx=SPlat_selectedSiteAndTx_all(trackID_selectedSiteAndTx_filtering == trackID_end);
        SPlon_selectedSiteAndTx=SPlon_selectedSiteAndTx_all(trackID_selectedSiteAndTx_filtering == trackID_end);
        SPAlt_selectedSiteAndTx=SPAlt_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        incAngle_selectedSiteAndTx=incAngle_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        PRN_selectedSiteAndTx=PRN_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        SVN_selectedSiteAndTx=SVN_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        trackID_selectedSiteAndTx=trackID_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        GNSS_Week_selectedSiteAndTx=GNSS_Week_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        GNSS_SoW_selectedSiteAndTx=GNSS_SoW_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        TxID_selectedSiteAndTx=TxID_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_x_selectedSiteAndTx=Tx_x_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_y_selectedSiteAndTx=Tx_y_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_z_selectedSiteAndTx=Tx_z_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_Vx_selectedSiteAndTx=Tx_Vx_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_Vy_selectedSiteAndTx=Tx_Vy_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Tx_Vz_selectedSiteAndTx=Tx_Vz_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_x_selectedSiteAndTx=Rx_x_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_y_selectedSiteAndTx=Rx_y_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_z_selectedSiteAndTx=Rx_z_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_Vx_selectedSiteAndTx=Rx_Vx_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_Vy_selectedSiteAndTx=Rx_Vy_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_Vz_selectedSiteAndTx=Rx_Vz_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        Rx_roll_selectedSiteAndTx=Rx_roll_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end); %Amirreza
        Rx_pitch_selectedSiteAndTx=Rx_pitch_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end); %Amirreza
        Rx_yaw_selectedSiteAndTx=Rx_yaw_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end); %Amirreza

        if OrbitWithReflectivity 
            ReflectivitydB_selectedSiteAndTx=ReflectivitydB_selectedSiteAndTx(trackID_selectedSiteAndTx_filtering == trackID_end);
        end

        %amir

    numberSelectedSP=abs(indexSP_end-indexSP_start)+1;
            
    for k=1:numberSelectedSP
         if indexSP_start<indexSP_end
             j=indexSP_start+k-1; %descending tracks
         else
             j=indexSP_end+k-1; %ascending tracks
         end
    

         SPlatlonAlt_selectedSP(k,1)=SPlat_selectedSiteAndTx(j);
         SPlatlonAlt_selectedSP(k,2)=SPlon_selectedSiteAndTx(j);
         SPlatlonAlt_selectedSP(k,3)=SPAlt_selectedSiteAndTx(j);
         incAngle_selectedSP(k)=incAngle_selectedSiteAndTx(j);
         Tx_selected(k,1)=Tx_x_selectedSiteAndTx(j);
         Tx_selected(k,2)=Tx_y_selectedSiteAndTx(j);
         Tx_selected(k,3)=Tx_z_selectedSiteAndTx(j);
         VTx_selected(k,1)=Tx_Vx_selectedSiteAndTx(j);
         VTx_selected(k,2)=Tx_Vy_selectedSiteAndTx(j);
         VTx_selected(k,3)=Tx_Vz_selectedSiteAndTx(j);
         Rx_selected(k,1)=Rx_x_selectedSiteAndTx(j);
         Rx_selected(k,2)=Rx_y_selectedSiteAndTx(j);
         Rx_selected(k,3)=Rx_z_selectedSiteAndTx(j);
         VRx_selected(k,1)=Rx_Vx_selectedSiteAndTx(j);
         VRx_selected(k,2)=Rx_Vy_selectedSiteAndTx(j);
         VRx_selected(k,3)=Rx_Vz_selectedSiteAndTx(j);
         GNSS_Week_selected(k)=GNSS_Week_selectedSiteAndTx(j);
         GNSS_SoW_selected(k)=GNSS_SoW_selectedSiteAndTx(j);
         trackID_selected(k)=trackID_selectedSiteAndTx(j); %%added by aneesha
         PRN_selected(k)=PRN_selectedSiteAndTx(j);%%added by aneesha
         SVN_selected(k)=SVN_selectedSiteAndTx(j);%%added by aneesha

         Rx_roll_selected(k)= Rx_roll_selectedSiteAndTx(j); % Added by amirreza
         Rx_pitch_selected(k)= Rx_pitch_selectedSiteAndTx(j); % Added by amirreza
         Rx_yaw_selected(k)= Rx_yaw_selectedSiteAndTx(j); % Added by amirreza

%          GnssBlock_selected(k)=GnssBlock_selectedSiteAndTx(j);
        if OrbitWithReflectivity
            ReflectivitydB_selected(k)=ReflectivitydB_selectedSiteAndTx(j); 
        end

         filenameFigDEMandDDM(k) = cellstr(['W' num2str(GNSS_Week_selected(k)) 'SoW' num2str(GNSS_SoW_selected(k)) 'T' num2str(trackID_selected(k))]);
    
    end

    if tag.plotSPsMap == 1
        %show land cover map in lat-lon
        %plot the selected SP over the land cover map
        figLandCoverAndSPs = figure('Name', 'Land cover map', 'NumberTitle','off','OuterPosition', [750 300 600 500],'Units','normalized');%[50 50 700 510]
        subplot('Position', [0.1 0.12 0.65 0.82]);
        mapshow(coverMap(:,:,2),coverMap(:,:,1),coverMap(:,:,3), 'DisplayType', 'texturemap')
        axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
        colormap(mapOfColorCover)
        colorbar('Ticks',[0,1,2,3,4,5,6],...
                'TickLabels',{'0 water bodies','1 broadleaved forest','2 niddleleaved forest','3 cropland/grassland/shrubland',...
                '4 bare areas','5 flooded forest','6 flooded shrubland'})
        caxis([0 6])
        title('Land cover map')
        xlabel('Longitude (deg)')
        ylabel('Latitude (deg)')
        mapshow(SPlatlonAlt_selectedSP(:,2),SPlatlonAlt_selectedSP(:,1),'DisplayType', 'point','Marker','o','MarkerFaceColor','w','MarkerEdgeColor','w')
        saveas(figLandCoverAndSPs, strjoin([outputDirectory 'landCoverAndSelectedSPs.fig'],''));

        %show the DEM
        %plot the selected SPs over the DEM
        if DEMinfo.flatsurface==0
            figDEMAndSPs = figure('Name', 'DEM', 'NumberTitle','off','OuterPosition', [150 300 600 500],'Units','normalized'); %added by ansha initial "figure"
            mapshow(DEM_lla(:,:,2),DEM_lla(:,:,1),DEM_lla(:,:,3), 'DisplayType', 'texturemap')
            axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
            cbarDEM=colorbar;
            cbarDEM.Label.String = 'Surface elevation (m)';
            title('DEM - SP white circle') 
            xlabel('Longitude (deg)')
            ylabel('Latitude (deg)')
            hold on
            scatter(SPlatlonAlt_selectedSP(:,2),SPlatlonAlt_selectedSP(:,1),25,[1,1,1],'filled')
            saveas(figDEMAndSPs, strjoin([outputDirectory 'DEMAndSelectedSPs.fig'],''));
        end
    end

else
%% Global mode or MonteCarlo mode (ALL TRACKS will be simulated)

    [length1,~]=size(trackID_selectedSiteAndTx);
    Num_SP = Sampling{3} ;
    if Num_SP > 0 && Num_SP<length1 

        indexSP_start=1;

        numberSelectedSP=abs(round(linspace(indexSP_start,length1,Num_SP)));

        SPlatlonAlt_selectedSP(:,1)=SPlat_selectedSiteAndTx(numberSelectedSP);
        SPlatlonAlt_selectedSP(:,2)=SPlon_selectedSiteAndTx(numberSelectedSP);
        SPlatlonAlt_selectedSP(:,3)=SPAlt_selectedSiteAndTx(numberSelectedSP);
        incAngle_selectedSP(:)=incAngle_selectedSiteAndTx(numberSelectedSP);
        Tx_selected(:,1)=Tx_x_selectedSiteAndTx(numberSelectedSP);
        Tx_selected(:,2)=Tx_y_selectedSiteAndTx(numberSelectedSP);
        Tx_selected(:,3)=Tx_z_selectedSiteAndTx(numberSelectedSP);
        VTx_selected(:,1)=Tx_Vx_selectedSiteAndTx(numberSelectedSP);
        VTx_selected(:,2)=Tx_Vy_selectedSiteAndTx(numberSelectedSP);
        VTx_selected(:,3)=Tx_Vz_selectedSiteAndTx(numberSelectedSP);
        Rx_selected(:,1)=Rx_x_selectedSiteAndTx(numberSelectedSP);
        Rx_selected(:,2)=Rx_y_selectedSiteAndTx(numberSelectedSP);
        Rx_selected(:,3)=Rx_z_selectedSiteAndTx(numberSelectedSP);
        VRx_selected(:,1)=Rx_Vx_selectedSiteAndTx(numberSelectedSP);
        VRx_selected(:,2)=Rx_Vy_selectedSiteAndTx(numberSelectedSP);
        VRx_selected(:,3)=Rx_Vz_selectedSiteAndTx(numberSelectedSP);
        GNSS_Week_selected(:)=GNSS_Week_selectedSiteAndTx(numberSelectedSP);
        GNSS_SoW_selected(:)=GNSS_SoW_selectedSiteAndTx(numberSelectedSP);
        trackID_selected(:)=trackID_selectedSiteAndTx(numberSelectedSP); %%added by aneesha
        PRN_selected(:)=PRN_selectedSiteAndTx(numberSelectedSP);%%added by aneesha
        SVN_selected(:)=SVN_selectedSiteAndTx(numberSelectedSP);%%added by aneesha
        
        Rx_roll_selected(:) = Rx_roll_selectedSiteAndTx(numberSelectedSP); % Added by Amirreza
        Rx_pitch_selected(:) = Rx_pitch_selectedSiteAndTx(numberSelectedSP); % Added by Amirreza
        Rx_yaw_selected(:) = Rx_yaw_selectedSiteAndTx(numberSelectedSP); % Added by Amirreza


%          GnssBlock_selected(:)=GnssBlock_selectedSiteAndTx(numberSelectedSP);
        if OrbitWithReflectivity
            ReflectivitydB_selected(:)=ReflectivitydB_selectedSiteAndTx(numberSelectedSP); 
        end

        for k=1:Num_SP
            filenameFigDEMandDDM(k) = cellstr(['W' num2str(GNSS_Week_selected(k)) 'SoW' num2str(GNSS_SoW_selected(k)) 'T' num2str(trackID_selected(k))]);    
        end
    else 
        if  Num_SP>length1
            Num_SP=length1;
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            uiwait(msgbox({'Warning: The requested number of points are more than those available.';...
                            'All available SPs will be simulated.'},CreateStruct))
        end

        SPlatlonAlt_selectedSP(:,1)=SPlat_selectedSiteAndTx;
        SPlatlonAlt_selectedSP(:,2)=SPlon_selectedSiteAndTx;
        SPlatlonAlt_selectedSP(:,3)=SPAlt_selectedSiteAndTx;
        incAngle_selectedSP=incAngle_selectedSiteAndTx;
        Tx_selected(:,1)=Tx_x_selectedSiteAndTx;
        Tx_selected(:,2)=Tx_y_selectedSiteAndTx;
        Tx_selected(:,3)=Tx_z_selectedSiteAndTx;
        VTx_selected(:,1)=Tx_Vx_selectedSiteAndTx;
        VTx_selected(:,2)=Tx_Vy_selectedSiteAndTx;
        VTx_selected(:,3)=Tx_Vz_selectedSiteAndTx;
        Rx_selected(:,1)=Rx_x_selectedSiteAndTx;
        Rx_selected(:,2)=Rx_y_selectedSiteAndTx;
        Rx_selected(:,3)=Rx_z_selectedSiteAndTx;
        VRx_selected(:,1)=Rx_Vx_selectedSiteAndTx;
        VRx_selected(:,2)=Rx_Vy_selectedSiteAndTx;
        VRx_selected(:,3)=Rx_Vz_selectedSiteAndTx;
        GNSS_Week_selected=GNSS_Week_selectedSiteAndTx;
        GNSS_SoW_selected=GNSS_SoW_selectedSiteAndTx;
        trackID_selected=trackID_selectedSiteAndTx;
        PRN_selected=PRN_selectedSiteAndTx;
        SVN_selected=SVN_selectedSiteAndTx;

        Rx_roll_selected = Rx_roll_selectedSiteAndTx; % Added by Amirreza
        Rx_pitch_selected = Rx_pitch_selectedSiteAndTx; % Added by Amirreza
        Rx_yaw_selected = Rx_yaw_selectedSiteAndTx; % Added by Amirreza


%     GnssBlock_selected=GnssBlock_selectedSiteAndTx;
        if OrbitWithReflectivity
            ReflectivitydB_selected=ReflectivitydB_selectedSiteAndTx; 
        end
        for k=1:length(SPlat_selectedSiteAndTx)
            filenameFigDEMandDDM(k) = cellstr(['W' num2str(GNSS_Week_selected(k)) 'SoW' num2str(GNSS_SoW_selected(k)) 'T' num2str(trackID_selected(k))]);
        end
    end

    if tag.plotSPsMap == 1
        if DEMinfo.flatsurface==0
            %open and read the DEM of the area of interest       
            [DEM_lla,~,~] = readDEMforPlot(DEMinfo.filenameDEM30arcsec,DEMinfo.SPcoodinatesLimits);
            %show DEM in lat-lon and overlap the SP tracks
           if tag.GUI == 0
            figDEMAndSPs = figure('Name', 'DEM', 'NumberTitle','off','OuterPosition', [150 300 600 500],'Units','normalized'); %[50 300 600 500]
            mapshow(DEM_lla(:,:,2),DEM_lla(:,:,1),DEM_lla(:,:,3), 'DisplayType', 'texturemap')
            axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
            cbarDEM=colorbar;
            cbarDEM.Label.String = 'Surface elevation (m)';
            title('DEM - SP white circle') 
            xlabel('Longitude (deg)')
            ylabel('Latitude (deg)')
            hold on
            scatter(SPlatlonAlt_selectedSP(:,2),SPlatlonAlt_selectedSP(:,1),25,[1,1,1],'filled')
            saveas(figDEMAndSPs,strjoin([outputDirectory 'DEMAndSelectedSPs.fig'],''));
           end
        end

        %open and read the land cover map of the area of interest
        SPlatlon_selectedSite(:,1)=SPlat_selectedSiteAndTx;
        SPlatlon_selectedSite(:,2)=SPlon_selectedSiteAndTx;
        [coverMap,~]=readCoverMapCCI(filenameCoverMap,SPlatlon_selectedSite,1.0);
        %plot the SP tracks over the land cover map of the area of interest
        mapOfColorCover=[0 0 1;0 1 0;1 0 1;0 1 1;1 0 0;1 1 0;0 0 0];
        if tag.GUI ==0
            figLandCoverAndSPs = figure('Name', 'Land cover map', 'NumberTitle','off','OuterPosition', [750 300 600 500],'Units','normalized'); %[900 300 600 500]
            subplot('Position', [0.1 0.12 0.65 0.82]);
            mapshow(coverMap(:,:,2),coverMap(:,:,1),coverMap(:,:,3), 'DisplayType', 'texturemap')
            axis([DEMinfo.SPcoodinatesLimits(3) DEMinfo.SPcoodinatesLimits(4) DEMinfo.SPcoodinatesLimits(1) DEMinfo.SPcoodinatesLimits(2)])
            colormap(mapOfColorCover)
            colorbar('Ticks',[0,1,2,3,4,5,6],...
                    'TickLabels',{'0 water bodies','1 broadleaved forest','2 needleleaved forest','3 cropland/grassland/shrubland',...
                    '4 bare areas','5 flooded forest','6 flooded shrubland'})
            caxis([0 6])
            title('Land cover map - SP white circle')
            xlabel('Longitude (deg)')
            ylabel('Latitude (deg)')
            hold on
            scatter(SPlatlonAlt_selectedSP(:,2),SPlatlonAlt_selectedSP(:,1),25,[1,1,1],'filled')
            saveas(figLandCoverAndSPs, strjoin([outputDirectory 'landCoverAndSelectedSPs.fig'],''));
        end
        if tag.GUI == 0
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            uiwait(msgbox({'Either Global mode or the Montecarlo mode has been activated.';...
                                'All SPs shown in the figure will be simulated.'},CreateStruct))
        end

    end
end

if ~OrbitWithReflectivity 
    ReflectivitydB_selected=NaN; 
end
if tag.GUI == 0
    ProceedQuestion  = questdlg(...
        ' Would you like to proceed with SAVERS simulations or to exit the program? ',...
        ' SP selection option ',...
        'Proceed with SAVERS simulations','Exit the program','Proceed with SAVERS simulations');
    if size(ProceedQuestion,1)<1
        errordlg('Exit before ending the simulation','Stop');
        Tx_selected=0; VTx_selected=0;Rx_selected=0;VRx_selected=0;SPlatlonAlt_selectedSP=[0,0,0];incAngle_selectedSP=0;
        PRN_selected=0;GNSS_Week_selected=0;GNSS_SoW_selected=0;trackID_selected=0;SVN_selected=0;goToSimulation=0;
        OrbitWithReflectivity=0;ReflectivitydB_selected=0;Rx_roll_selected=0;Rx_pitch_selected=0;Rx_yaw_selected=0;
        filenameFigDEMandDDM='';
        return
    end
    
    switch ProceedQuestion
        case 'Proceed with SAVERS simulations'
            goToSimulation=1;
            close all
        case 'Exit the program'
            goToSimulation=0;
            disp('.......................................')
            disp(' ### SIMULATION STOPPED BY THE USER ### ')
            disp('........................................')
            close all
            return
    end
else 
        goToSimulation=1;
        close all
end