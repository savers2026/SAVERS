%% ************************************************************************
% MODULE NAME:      write_L1a_metadata.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      Aneesha Musunuri, N. Pierdicca
% @copyright:   Sapienza University of Rome
% Original version:      Sep 2021 - Sep 2022 by Aneesha Musunuri, N. Pierdicca
% Main updates:          Dec 2022 by Laura Dente and Davide Comite
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module writes the metadata_1a.nc L1a product of the Hydrognss
% E2ES.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function write_L1a_metadata(satellite,TRAp,GNSS_SoW,GNSS_Week,Rx,VRx,Tx,VTx,PRN,TrackID,SVN,name_data_tag,geoSYSp,...
                                IntegrationMidPointTime,ReflectionHeight,SatSpecularPointLat,SatSpecularPointLon,...
                                OnBoardDelayReferencedToDirect_1,OnBoardDelayReferencedToDirect_5,...
                                max_sigPower_Noise_LR1,max_sigPower_Noise_LR5,max_sigPower_Noise_RR1,max_sigPower_Noise_RR5,...
                                outpathname,foldername, IntegrationMidPointTime_coe)

[v, w] = unique( TrackID, 'stable' );
totaloftrackid=numel(TrackID);
nooftracks = length(w);
j = w(2:nooftracks,1)-1;
k=numel(j);
for l=1:nooftracks
    if l <= k
        sizeoftracks(l)=w(l+1)-w(l);
    else
        sizeoftracks(l)=(totaloftrackid+1)-w(l);
    end
end
newtrackid=v;
%global variables from orbitdata
% ReceiverMetaDataPerEpoch='N_epoch';
% IntegrationStartTimeGPSWeeks=GNSS_Week;
% IntegrationStartTimeGPSSeconds=GNSS_SoW;
IntegrationMidTimeGPSWeeks=GNSS_Week;
IntegrationMidTimeGPSSeconds=GNSS_SoW;
ReceiverPositionX=Rx(:,1);
ReceiverPositionY=Rx(:,2);
ReceiverPositionZ=Rx(:,3);
ReceiverVelocityX=VRx(:,1);
ReceiverVelocityY=VRx(:,2);
ReceiverVelocityZ=VRx(:,3);
GNSSBroadcastPositionX=Tx(:,1);
GNSSBroadcastPositionY=Tx(:,2);
GNSSBroadcastPositionZ=Tx(:,3);
GNSSBroadcastVelocityX=VTx(:,1);
GNSSBroadcastVelocityY=VTx(:,2);
GNSSBroadcastVelocityZ=VTx(:,3);
InEclipse=zeros(size(ReceiverVelocityZ));

AttitudeUncertainty_roll   = repmat(TRAp.RxAttitude_angles_GUI(1), size(ReceiverVelocityZ));
AttitudeUncertainty_pitch  = repmat(TRAp.RxAttitude_angles_GUI(2), size(ReceiverVelocityZ));
AttitudeUncertainty_yaw    = repmat(TRAp.RxAttitude_angles_GUI(3), size(ReceiverVelocityZ));

Attitude_roll              = repmat(TRAp.Rx_nomAtt_angles(1),size(ReceiverVelocityZ));
Attitude_pitch             = repmat(TRAp.Rx_nomAtt_angles(2),size(ReceiverVelocityZ));
Attitude_yaw               = repmat(TRAp.Rx_nomAtt_angles(3),size(ReceiverVelocityZ));

%global variables from metadatafile
% IntegrationMidPointTime = datenum(1980, 1, 6, 0 , 0 ,0)+GNSS_Week.*7+GNSS_SoW./86400;
% IntegrationStartTime = ones(size(ReceiverVelocityZ));%dummyfornow
% IntegrationEndTime   = ones(size(ReceiverVelocityZ));%dummyfornow
IntegrationStartTime = IntegrationMidPointTime - geoSYSp.incoIntTime_sec/2;
IntegrationEndTime   = IntegrationMidPointTime + geoSYSp.incoIntTime_sec/2;

ClockStatus    = uint8(zeros(size(ReceiverVelocityZ)));%uncorrected
ClockErrorRate = repmat(1574.42,size(ReceiverVelocityZ));

GDOP = repmat(2,size(ReceiverVelocityZ));

AttitudeSource      = int16(ones(size(ReceiverVelocityZ))); %platform
InternalTemperature = geoSYSp.TRx_K - 273.15;

[InternalTemperatures0,InternalTemperatures1,InternalTemperatures2,InternalTemperatures3,InternalTemperatures4,InternalTemperatures5,InternalTemperatures6] = deal(repmat(InternalTemperature,(size(ReceiverVelocityZ))));

ActiveTracks0 = ones(size(ReceiverVelocityZ));
ActiveTracks1 = ones(size(ReceiverVelocityZ));
ActiveTracks2 = ones(size(ReceiverVelocityZ));
ActiveTracks3 = ones(size(ReceiverVelocityZ));

TimeDiffToClosestPVT    = zeros(size(ReceiverVelocityZ));
ReceiverPositionSource  = uint8(zeros(size(ReceiverVelocityZ)));%predicted
StarTrackerBlindingFlag = zeros(size(ReceiverVelocityZ)); %dummy
NumberOfTracks          = repmat(4,size(ReceiverVelocityZ));
NumberOfPolarisations   = repmat(2,size(ReceiverVelocityZ));

if max_sigPower_Noise_LR5(1) == 0
    NumberOfFrequencies         = ones(size(ReceiverVelocityZ));
    noofchannels                = 2;
    ReflectionChannelNumberList = uint16([1,2]);
else
    NumberOfFrequencies         = ones(size(ReceiverVelocityZ))+1;
    noofchannels                = 4;
    ReflectionChannelNumberList = uint16([1,2,3,4]);
end

%% lla to xyz
wgs84   = wgs84Ellipsoid('meters');
[X,Y,Z] = geodetic2ecef(wgs84,SatSpecularPointLat,SatSpecularPointLon,ReflectionHeight);

[lats,lons,alts] = ecef2geodetic(wgs84,ReceiverPositionX,ReceiverPositionY,ReceiverPositionZ);

ReceiverSubSatLatitude    = lats;
ReceiverSubSatLongitude   = lons;
SatSpecularPointPositionX = X;
SatSpecularPointPositionY = Y;
SatSpecularPointPositionZ = Z;
track_id                  = uint32(TrackID);

PRN_1 = PRN;
SVN_1 = SVN;

DirectSignalInDDM=zeros(size(ReceiverVelocityZ));

[ADCPercentagePositive,ADCOffset,ADCPercentageMagnitude] = deal(zeros(size(ReceiverVelocityZ)));
ncid         = netcdf.create(strjoin([outpathname '\' 'metadata_1a.nc'],''),'netcdf4');

% FirstTimeStamp=datetime('now','Format','dd-MM-uuuu HH:mm:ss');
% FirstTimeStamp=datestr(FirstTimeStamp);
% LastTimeStamp=datetime('now','Format','dd-MM-uuuu HH:mm:ss');
% LastTimeStamp=datestr(LastTimeStamp);
FirstTimeStamp = datestr(datetime(IntegrationMidPointTime(1), 'ConvertFrom', 'datenum','Format','dd-MM-uuuu HH:mm:ss'));
LastTimeStamp  = datestr(datetime(IntegrationMidPointTime(end), 'ConvertFrom', 'datenum','Format','dd-MM-uuuu HH:mm:ss'));
%dimensions
% IntegrationMidPointTime=IntegrationMidPointTime';
[a b] = size(IntegrationMidPointTime);
dimid = netcdf.defDim(ncid,'DateTimes',a); %4.6hours (timeperiod)
ProcessingType = "Ground4";
if satellite == "HydroGNSS-1"
    SourceNumber   = uint32(1);
else
    SourceNumber   = uint32(2);
end
L1a_ProcessorVersion    = string(1.00);
FPGAVersion             = single(0.20);
SoftwareReceiverVersion = single(2.0);

IncoherentMeasurementFrequency = single(1/geoSYSp.incoIntTime);

GpsUTCOffset = int16(18);
PtrId        = int32(0);
MasterNavChannelNumber      = uint16(4);
GnssBlock                   = int32(7);
GnssBlock_units             = 'GnssBlock_units';
AllocationMode              = int32(0);
AllocationMode_units        = 'AllocationMode_units';
%SignalFrequency             = geoSYSp.TxSignalButton;
% GNSSConstellation=geoSYSp.indexTxSignal;

if geoSYSp.indexTxSignal == 0 || geoSYSp.indexTxSignal ==1
    GNSSConstellation = uint32(0);
    SignalFrequency   = 'GPS_L';
else
    GNSSConstellation = uint32(2);
    SignalFrequency   = 'Galileo_E';
end

% SignalPolarisation='LHCP';
GNSSConstellation_units = '{0:GPS,2:GALILEO}';
% FrontEnd                = int32(1);
FrontEnd_units          = '{ 0:Nadir_L1E1_LHCP, 1:Nadir_L1E1_RHCP, 2:Nadir_L5E5_LHCP, 3:Nadir_L5E5_RHCP, 4:ZenithL_L1E1_RHCP,5:ZenithL_L1E1_LHCP(N/A), 6:ZenithL_L5E5_RHCP, 7:ZenithL_L5E5_LHCP(N/A)}';
SamplingFrequency       = double(geoSYSp.SamplingFrequency_Hz);
TrackingMode            = int32(1);
TrackingMode_units      = '{1:DirectSignalOverride, 2:OpenLoopReflection }';
FrontEndGainMode        = int32(0);
FrontEndGainMode_units  = '{0:FixedGain, 1:AutoUnmonitored, 2:AutomaticMonitored}';
FrontEndGainControlSetting = uint32(36);
LNASwitchStatus            = int32(0);
LNASwitchStatus_units      = '{0:Antenna,1:Reference}';
NumberOfDelayPixels        = uint32(geoSYSp.nx_Rdelay);
NumberOfDopplerPixels      = uint32(geoSYSp.ny_Dfreq);
DelayResolution            = double(10^3.*geoSYSp.mpas); 
CodeDelaySpacingSamplesBetweenPixels = uint16(geoSYSp.CodeDelaySpacingSamplesBetweenPixels);
DopplerResolution                    = double(geoSYSp.fpas);
IncoherentIntegrations               = uint32(geoSYSp.incoIntTime);
%TrackingOffsetDelayInPixels          = double(floor(delayOffset/DelayResolution*1000));
TrackingOffsetDelayInPixels          = geoSYSp.delayOffset_bin;

% TrackingOffsetDelayNs=double(15*(geoSYSp.Dly_RNGmax/(NumberOfDelayPixels)*10^6));
% DelayResolution=double(double(CodeDelaySpacingSamplesBetweenPixels)/SamplingFrequency*1e9);
%TrackingOffsetDelayNs         = double(TrackingOffsetDelayInPixels* DelayResolution);
TrackingOffsetDelayNs         = geoSYSp.delayOffset_microsec*10^3;
TrackingOffsetDopplerHz       = double(0.0);
TrackingOffsetDopplerInPixels = double(0.0);
% CoherentIntegrationTime=uint32(geoSYSp.cohIntTime);
CoherentIntegrationTime       = double(geoSYSp.cohIntTime);
datecreated                   = datetime('now','Format','yyyy-MM-dd''T''hh:mm:ss');

%global attributes
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'Name','HydroGNSS L1a Metadata File');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'date_created',string(datecreated));
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'License','HydroGNSS GNSS-R Dataset');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'FileIDCode',foldername);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'DataTag',name_data_tag);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceConstellation','HYDROGNSS');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceNumber',SourceNumber);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'ProcessingType',ProcessingType);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'L1a_ProcessorVersion',L1a_ProcessorVersion);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'FPGAVersion',FPGAVersion);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SoftwareReceiverVersion',SoftwareReceiverVersion);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'IncoherentMeasurementFrequency',IncoherentMeasurementFrequency);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'GpsUTCOffset',GpsUTCOffset);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'FirstTimeStamp',FirstTimeStamp);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'LastTimeStamp',LastTimeStamp);

%global variables
% dimid1=netcdf.defDim(ncid,'x1',1);
% var1=netcdf.defVar(ncid,'ReceiverMetaDataPerEpoch','NC_STRING',dimid1);
% netcdf.putVar(ncid,var1,ReceiverMetaDataPerEpoch);
var2=netcdf.defVar(ncid,'IntegrationMidPointTime','NC_DOUBLE',dimid);%(gps_week +gps_sow ) from previous metadatafile
netcdf.putVar(ncid,var2,IntegrationMidPointTime(1:a));
var3=netcdf.defVar(ncid,'IntegrationMidPointTimeGPSWeeks','NC_SHORT',dimid);%gps_week in seconds (to recall - first week to be checked)
netcdf.putVar(ncid,var3,IntegrationMidTimeGPSWeeks(1:a));
var4=netcdf.defVar(ncid,'IntegrationMidPointTimeGPSSeconds','NC_DOUBLE',dimid);%gps_sow-0.5(incoherentintegrationtime)
netcdf.putVar(ncid,var4,IntegrationMidTimeGPSSeconds(1:a));
var5=netcdf.defVar(ncid,'IntegrationStartTime','NC_DOUBLE',dimid);%(gps_week +gps_sow ) - 0.5( )
netcdf.putVar(ncid,var5,IntegrationStartTime(1:a));
var6=netcdf.defVar(ncid,'IntegrationEndTime','NC_DOUBLE',dimid);%(gps_week +gps_sow ) + 0.5(incoherentintegrationtime)
netcdf.putVar(ncid,var6,IntegrationEndTime(1:a));
var7=netcdf.defVar(ncid,'ClockStatus','NC_UBYTE',dimid);
netcdf.putVar(ncid,var7,ClockStatus(1:a));
var8=netcdf.defVar(ncid,'ClockErrorRate','NC_FLOAT',dimid);
netcdf.putVar(ncid,var8,ClockErrorRate(1:a));
var9_1=netcdf.defVar(ncid,'ReceiverPositionX','NC_DOUBLE',dimid); %pos_vel_rx (1-3)
netcdf.putVar(ncid,var9_1,ReceiverPositionX(1:a));
var9_2=netcdf.defVar(ncid,'ReceiverPositionY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var9_2,ReceiverPositionY(1:a));
var9_3=netcdf.defVar(ncid,'ReceiverPositionZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var9_3,ReceiverPositionZ(1:a));
var10_1=netcdf.defVar(ncid,'ReceiverVelocityX','NC_DOUBLE',dimid); %pos_vel_rx(4-6)
netcdf.putVar(ncid,var10_1,ReceiverVelocityX(1:a));
var10_2=netcdf.defVar(ncid,'ReceiverVelocityY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var10_2,ReceiverVelocityY(1:a));
var10_3=netcdf.defVar(ncid,'ReceiverVelocityZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var10_3,ReceiverVelocityZ(1:a));
var11=netcdf.defVar(ncid,'ReceiverPositionSource','NC_UBYTE',dimid);
netcdf.putVar(ncid,var11,ReceiverPositionSource(1:a));
var11_1=netcdf.defVar(ncid,'ReceiverSubSatLatitude','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var11_1,ReceiverSubSatLatitude(1:a));
var11_2=netcdf.defVar(ncid,'ReceiverSubSatLongitude','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var11_2,ReceiverSubSatLongitude(1:a));
var12=netcdf.defVar(ncid,'GDOP','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var12,GDOP(1:a));
var13_1=netcdf.defVar(ncid,'AttitudeRoll','NC_FLOAT',dimid);% taken from old metadata file(zeros)
netcdf.putVar(ncid,var13_1,Attitude_roll(1:a));
var13_2=netcdf.defVar(ncid,'AttitudePitch','NC_FLOAT',dimid);%added
netcdf.putVar(ncid,var13_2,Attitude_pitch(1:a));
var13_3=netcdf.defVar(ncid,'AttitudeYaw','NC_FLOAT',dimid);%added
netcdf.putVar(ncid,var13_3,Attitude_yaw(1:a));
var14_1=netcdf.defVar(ncid,'AttitudeUncertaintyRoll','NC_FLOAT',dimid);% taken from old metadata file(zeros)
netcdf.putVar(ncid,var14_1,AttitudeUncertainty_roll(1:a));
var14_2=netcdf.defVar(ncid,'AttitudeUncertaintyPitch','NC_FLOAT',dimid);
netcdf.putVar(ncid,var14_2,AttitudeUncertainty_pitch(1:a));
var14_3=netcdf.defVar(ncid,'AttitudeUncertaintyYaw','NC_FLOAT',dimid);
netcdf.putVar(ncid,var14_3,AttitudeUncertainty_yaw(1:a));
var15=netcdf.defVar(ncid,'StarTrackerBlindingFlag','NC_USHORT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var15,StarTrackerBlindingFlag(1:a));
var16=netcdf.defVar(ncid,'AttitudeSource','NC_SHORT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var16,AttitudeSource(1:a));
var17=netcdf.defVar(ncid,'InEclipse','NC_SHORT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var17,InEclipse(1:a));
var18=netcdf.defVar(ncid,'InternalTemperatures0','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18,InternalTemperatures0(1:a));
var18_a=netcdf.defVar(ncid,'InternalTemperatures1','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_a,InternalTemperatures1(1:a));
var18_b=netcdf.defVar(ncid,'InternalTemperatures2','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_b,InternalTemperatures2(1:a));
var18_c=netcdf.defVar(ncid,'InternalTemperatures3','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_c,InternalTemperatures3(1:a));
var18_d=netcdf.defVar(ncid,'InternalTemperatures4','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_d,InternalTemperatures4(1:a));
var18_e=netcdf.defVar(ncid,'InternalTemperatures5','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_e,InternalTemperatures5(1:a));
var18_f=netcdf.defVar(ncid,'InternalTemperatures6','NC_FLOAT',dimid);% taken from old metadata file
netcdf.putVar(ncid,var18_f,InternalTemperatures6(1:a));
% var19=netcdf.defVar(ncid,'PlatformTemperatures','NC_DOUBLE',dimid);% taken from old metadata file
% netcdf.putVar(ncid,var19,PlatformTemperatures(1:a));
var20=netcdf.defVar(ncid,'ActiveTracks0','NC_INT',dimid);% 4 tracks for examples
netcdf.putVar(ncid,var20,ActiveTracks0(1:a));
var20_a=netcdf.defVar(ncid,'ActiveTracks1','NC_INT',dimid);% 4 tracks for examples
netcdf.putVar(ncid,var20_a,ActiveTracks1(1:a));
var20_b=netcdf.defVar(ncid,'ActiveTracks2','NC_INT',dimid);% 4 tracks for examples
netcdf.putVar(ncid,var20_b,ActiveTracks2(1:a));
var20_c=netcdf.defVar(ncid,'ActiveTracks3','NC_INT',dimid);% 4 tracks for examples
netcdf.putVar(ncid,var20_c,ActiveTracks3(1:a));
var20_d=netcdf.defVar(ncid,'TimeDiffToClosestPVT','NC_FLOAT',dimid);
netcdf.putVar(ncid,var20_d,TimeDiffToClosestPVT(1:a));
var21=netcdf.defVar(ncid,'NumberOfTracks','NC_SHORT',dimid);
netcdf.putVar(ncid,var21,NumberOfTracks(1:a));
var22=netcdf.defVar(ncid,'NumberOfPolarisations','NC_SHORT',dimid);
netcdf.putVar(ncid,var22,NumberOfPolarisations(1:a));
var23=netcdf.defVar(ncid,'NumberOfFrequencies','NC_SHORT',dimid);
netcdf.putVar(ncid,var23,NumberOfFrequencies(1:a));
% 
for  n=1:nooftracks % no of trackgroups
    if n<=10
        trackid(n) = netcdf.defGrp(ncid,['00000' num2str(n-1)]);
    elseif n>=11 && n<=100 
        trackid(n) = netcdf.defGrp(ncid,['0000' num2str(n-1)]);
    elseif n>=101 && n<=1000
        trackid(n) = netcdf.defGrp(ncid,['000' num2str(n-1)]);
    elseif n>=1001 && n<=10000
        trackid(n) = netcdf.defGrp(ncid,['00' num2str(n-1)]);
    elseif n>=10001 && n<=100000
        trackid(n) = netcdf.defGrp(ncid,['0' num2str(n-1)]);
    else
        trackid(n) = netcdf.defGrp(ncid,num2str(n-1));
    end
    %index=find(n==track_id(:,1)|n==track_id(:,2)|n==track_id(:,3)|n==track_id(:,4));
    %[a1 b1]=size(index);
    index = find(newtrackid(n)==TrackID);
    [l1 l2]=size(index);
    dimid3 = netcdf.defDim(trackid(n),'DateTimeTrack',l1);
    % modifica Mauro. Inserito Uint32((n-1)) al posto di newtrackid(n) 
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'TrackID',uint32((n-1)));
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'TrackIDOrbit',uint32(track_id(index(1))));
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'DataTag','data_tag');
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'PtrId',PtrId);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'ReflectionChannelNumberList',ReflectionChannelNumberList);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'MasterNavChannelNumber',MasterNavChannelNumber);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'GNSSConstellation',GNSSConstellation);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'GNSSConstellation_units',GNSSConstellation_units);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'PRN',uint16((PRN_1(index(1)))));
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'SVN',uint16((SVN_1(index(1)))));
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'GnssBlock',GnssBlock);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'GnssBlock_units',GnssBlock_units);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'AllocationMode',AllocationMode);
    netcdf.putAtt(trackid(n),netcdf.getConstant("NC_GLOBAL"),'AllocationMode_units',AllocationMode_units);
    var1_tr = netcdf.defVar(trackid(n),'IntegrationMidPointTime','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),var1_tr,IntegrationMidPointTime(index));
    dimid_ln=netcdf.defDim(trackid(n),'Temp Axis',4);
    dimid_n=([dimid_ln dimid3]);
    var2_tr1=netcdf.defVar(trackid(n),'Temp Axis','NC_DOUBLE',dimid_ln);
    netcdf.putVar(trackid(n),var2_tr1,1:4);
    %           var2_tr = netcdf.defVar(trackid(n),'LNATemperatures','NC_FLOAT',dimid_n);
    %           netcdf.putVar(trackid(n),var2_tr,LNATemperatures(index,:));
    var3_tr = netcdf.defVar(trackid(n),'OnBoardReflectionHeight','NC_FLOAT',dimid3); % changed from NC_DOUBLE
    netcdf.putVar(trackid(n),var3_tr,ReflectionHeight(index));
    vara_tr=netcdf.defVar(trackid(n),'GNSSBroadcastPositionX','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),vara_tr,GNSSBroadcastPositionX(index));
    varb_tr=netcdf.defVar(trackid(n),'GNSSBroadcastPositionY','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),varb_tr,GNSSBroadcastPositionY(index));
    varc_tr=netcdf.defVar(trackid(n),'GNSSBroadcastPositionZ','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),varc_tr,GNSSBroadcastPositionZ(index));
    vara1_tr=netcdf.defVar(trackid(n),'GNSSBroadcastVelocityX','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),vara1_tr,GNSSBroadcastVelocityX(index));
    varb1_tr=netcdf.defVar(trackid(n),'GNSSBroadcastVelocityY','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),varb1_tr,GNSSBroadcastVelocityY(index));
    varc1_tr=netcdf.defVar(trackid(n),'GNSSBroadcastVelocityZ','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),varc1_tr,GNSSBroadcastVelocityZ(index));
    var4_tr = netcdf.defVar(trackid(n),'OnBoardSpecularPointPositionX','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),var4_tr,SatSpecularPointPositionX(index));
    var5_tr = netcdf.defVar(trackid(n),'OnBoardSpecularPointPositionY','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),var5_tr,SatSpecularPointPositionY(index));
    var6_tr = netcdf.defVar(trackid(n),'OnBoardSpecularPointPositionZ','NC_DOUBLE',dimid3);
    netcdf.putVar(trackid(n),var6_tr,SatSpecularPointPositionZ(index));
    var7_tr= netcdf.defVar(trackid(n),'OnBoardSpecularPointLat','NC_FLOAT',dimid3);
    netcdf.putVar(trackid(n),var7_tr,SatSpecularPointLat(index));
    var8_tr = netcdf.defVar(trackid(n),'OnBoardSpecularPointLon','NC_FLOAT',dimid3);
    netcdf.putVar(trackid(n),var8_tr,SatSpecularPointLon(index));
    var9_tr = netcdf.defVar(trackid(n),'OnBoardDelayReferencedToDirect','NC_FLOAT',dimid3);
    netcdf.putVar(trackid(n),var9_tr,OnBoardDelayReferencedToDirect_1(index));

    for m=1:noofchannels % no of channels
        trackid_ch(m) = netcdf.defGrp(trackid(n),['Channel' num2str(m)]);
        dimidch = netcdf.defDim(trackid_ch(m),'a1ch',l1);
        var9_ch = netcdf.defVar(trackid_ch(m),'IntegrationMidPointTime','NC_DOUBLE',dimidch);
        netcdf.putVar(trackid_ch(m),var9_ch,IntegrationMidPointTime(index));
        var10_ch = netcdf.defVar(trackid_ch(m),'DirectSignalInDDM','NC_UBYTE',dimidch);
        netcdf.putVar(trackid_ch(m),var10_ch,DirectSignalInDDM(index));
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'ChannelID',uint16(m));
        netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'TrackID',uint32((n-1)));
%        netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'SignalFrequency',replace(SignalFrequency,' ','_'));
        if m==1 || m==2
            netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'SignalFrequency',[SignalFrequency '1']);
        else
            netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'SignalFrequency',[SignalFrequency '5']);
        end
        if mod(m,2)~=0
            SignalPolarisation='LHCP';
        else
            SignalPolarisation='RHCP';
        end
        if m==1 || m==2 && SignalPolarisation == "LHCP"
            FrontEnd = int32(0);
        elseif m==1 || m==2 && SignalPolarisation =="RHCP"
            FrontEnd = int32(1);
        elseif m>2 && SignalPolarisation =="LHCP"   
             FrontEnd = int32(2);
        elseif m>2 && SignalPolarisation =="RHCP"
            FrontEnd = int32(3);
        end

        netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'SignalPolarisation',SignalPolarisation);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'FrontEnd',FrontEnd);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'FrontEnd_units',FrontEnd_units);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'SamplingFrequency',SamplingFrequency);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'TrackingMode',TrackingMode);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'TrackingMode_units',TrackingMode_units);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'FrontEndGainMode',FrontEndGainMode);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'FrontEndGainMode_units',FrontEndGainMode_units);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'FrontEndGainControlSetting',FrontEndGainControlSetting);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'LNASwitchStatus',LNASwitchStatus);
        netcdf.putAtt(trackid_ch(m),netcdf.getConstant("NC_GLOBAL"),'LNASwitchStatus_units',LNASwitchStatus_units);

        trackid_n1(m)= netcdf.defGrp(trackid_ch(m),'Incoherent');
        dimidchn1 = netcdf.defDim(trackid_n1(m),'n1',l1);
        if m==1 || m==2
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'DelayResolution',DelayResolution(1));
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'CodeDelaySpacingSamplesBetweenPixels',CodeDelaySpacingSamplesBetweenPixels(1));
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'TrackingOffsetDelayNs',TrackingOffsetDelayNs(1));
        else
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'DelayResolution',DelayResolution(2));
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'CodeDelaySpacingSamplesBetweenPixels',CodeDelaySpacingSamplesBetweenPixels(2));
            netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'TrackingOffsetDelayNs',TrackingOffsetDelayNs(2));
        end
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'DopplerResolution',DopplerResolution);
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'IncoherentIntegrations',IncoherentIntegrations);
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'NumberOfDelayPixels',NumberOfDelayPixels);
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'NumberOfDopplerPixels',NumberOfDopplerPixels);        
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'TrackingOffsetDelayInPixels',TrackingOffsetDelayInPixels);
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'TrackingOffsetDopplerHz',TrackingOffsetDopplerHz);
        netcdf.putAtt(trackid_n1(m),netcdf.getConstant("NC_GLOBAL"),'TrackingOffsetDopplerInPixels',TrackingOffsetDopplerInPixels);
        var11_n1=netcdf.defVar(trackid_n1(m),'IntegrationMidPointTime','NC_DOUBLE',dimidchn1);
        netcdf.putVar(trackid_n1(m),var11_n1,IntegrationMidPointTime(index));
        var12_n1=netcdf.defVar(trackid_n1(m),'ADCPercentagePositive','NC_FLOAT',dimidchn1);
        netcdf.putVar(trackid_n1(m),var12_n1,ADCPercentagePositive(index));
        var13_n1=netcdf.defVar(trackid_n1(m),'ADCOffset','NC_FLOAT',dimidchn1);
        netcdf.putVar(trackid_n1(m),var13_n1,ADCOffset(index));
        var14_n1=netcdf.defVar(trackid_n1(m),'ADCPercentageMagnitude','NC_FLOAT',dimidchn1);
        netcdf.putVar(trackid_n1(m),var14_n1,ADCPercentageMagnitude(index))
        var16_n1=netcdf.defVar(trackid_n1(m),'DDMOutputNumericalScaling','NC_DOUBLE',dimidchn1);
        switch m
            case 1
                DDMOutputNumericalScaling = max_sigPower_Noise_LR1;
            case 2
                DDMOutputNumericalScaling = max_sigPower_Noise_RR1;
            case 3
                DDMOutputNumericalScaling = max_sigPower_Noise_LR5;
            case 4
                DDMOutputNumericalScaling = max_sigPower_Noise_RR5;            
        end
        netcdf.putVar(trackid_n1(m),var16_n1,DDMOutputNumericalScaling(index));

        trackid_n2(m)= netcdf.defGrp(trackid_ch(m),'Coherent');
        %modified by Laura
%         dimidchn2 = netcdf.defDim(trackid_n2(m),'n12',l1);
%         var15_n2=netcdf.defVar(trackid_n2(m),'IntegrationMidPointTime','NC_DOUBLE',dimidchn2);
%         netcdf.putVar(trackid_n2(m),var15_n2,IntegrationMidPointTime(index));
%         netcdf.putAtt(trackid_n2(m),netcdf.getConstant("NC_GLOBAL"),'CoherentIntegrationTime',CoherentIntegrationTime);

        dimidchn2 = netcdf.defDim(trackid_n2(m),'n2',length(IntegrationMidPointTime_coe));
        var15_n2=netcdf.defVar(trackid_n2(m),'IntegrationMidPointTime','NC_DOUBLE',dimidchn2);
        netcdf.putVar(trackid_n2(m),var15_n2,IntegrationMidPointTime_coe);
        netcdf.putAtt(trackid_n2(m),netcdf.getConstant("NC_GLOBAL"),'CoherentIntegrationTime',CoherentIntegrationTime);
    end

end

netcdf.close(ncid);
% Metadata_HydroGNSS =ncinfo('metadata_1a.nc');
end

