%% ************************************************************************
% MODULE NAME:      write_L1a_directsignalpower.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      Aneesha Musunuri, N. Pierdicca
% @copyright:   Sapienza University of Rome
% Original version:      Sep 2021 - Sep 2022 by Aneesha Musunuri, N. Pierdicca
% Main updates:          Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module writes the directSignalPower.nc L1a product of the Hydrognss
% E2ES.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function write_L1a_DSpow(TRAp,Gsys_dB,Rx,VRx,Tx,VTx,SignalPower_ds,NoisePower_ds,BlackBodyPower_ds,PRN,SVN,IntegrationMidPointTime,geoSYSp,outpathname)

SourceNumber=uint32(1);
L1a_ProcessorVersion=single(1.00);
[a,~]=size(IntegrationMidPointTime);
StartTime=datetime(IntegrationMidPointTime(1), 'ConvertFrom', 'datenum');
StopTime=datetime(IntegrationMidPointTime(end), 'ConvertFrom', 'datenum');
% DateTime=string(datetime(IntegrationMidPointTime,'ConvertFrom','datenum'));
DateTime=IntegrationMidPointTime;
GnssBlock=repmat(7,size(IntegrationMidPointTime));
FrontEndGainMode=zeros(size(IntegrationMidPointTime));
AzimuthReceiverAntenna=zeros(size(IntegrationMidPointTime));
ElevationReceiverAntenna=zeros(size(IntegrationMidPointTime));
ElevationTransmitterAntenna=zeros(size(IntegrationMidPointTime));
SignalPower=SignalPower_ds;
NoisePower=NoisePower_ds;
BlackBodyPower=BlackBodyPower_ds;
SystemGainBB=repmat(Gsys_dB,size(IntegrationMidPointTime));
SystemGainBBComp=zeros(size(IntegrationMidPointTime));
RxTemperature=repmat(290,size(IntegrationMidPointTime));
InEclipse=zeros(size(IntegrationMidPointTime));
ReceiverPositionX=Rx(:,1);
ReceiverPositionY=Rx(:,2);
ReceiverPositionZ=Rx(:,3);
ReceiverVelocityX=VRx(:,1);
ReceiverVelocityY=VRx(:,2);
ReceiverVelocityZ=VRx(:,3);
TransmitterPositionX=Tx(:,1);
TransmitterPositionY=Tx(:,2);
TransmitterPositionZ=Tx(:,3);
TransmitterVelocityX=VTx(:,1);
TransmitterVelocityY=VTx(:,2);
TransmitterVelocityZ=VTx(:,3);
AttitudeRoll=repmat(TRAp.RxEuler_angles(1,1),size(ReceiverVelocityZ));
AttitudePitch=repmat(TRAp.RxEuler_angles(1,2),size(ReceiverVelocityZ));
AttitudeYaw=repmat(TRAp.RxEuler_angles(1,3),size(ReceiverVelocityZ));
InternalTemperature=geoSYSp.TRx_K - 273.15;
[InternalTemperatures0,InternalTemperatures1,InternalTemperatures2,InternalTemperatures3,InternalTemperatures4,InternalTemperatures5,InternalTemperatures6]=deal(repmat(InternalTemperature,(size(ReceiverVelocityZ))));

SPAzimuthARF=zeros(size(IntegrationMidPointTime));
SPElevationARF=zeros(size(IntegrationMidPointTime));


ncid = netcdf.create(strjoin([outpathname '\' 'directSignalPowerL1E1.nc'],''),'netcdf4');
dimid = netcdf.defDim(ncid,'DateTime',a);


var1=netcdf.defVar(ncid,'DateTime','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var1,DateTime);
var2=netcdf.defVar(ncid,'PRN','NC_USHORT',dimid);
netcdf.putVar(ncid,var2,uint16(PRN));
var3=netcdf.defVar(ncid,'SVN','NC_USHORT',dimid);
netcdf.putVar(ncid,var3,uint16(SVN));
var4=netcdf.defVar(ncid,'GnssBlock','NC_INT',dimid);
netcdf.putVar(ncid,var4,GnssBlock);
var5=netcdf.defVar(ncid,'FrontEndGainMode','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var5,FrontEndGainMode);
var6=netcdf.defVar(ncid,'AzimuthReceiverAntenna','NC_FLOAT',dimid);
netcdf.putVar(ncid,var6,AzimuthReceiverAntenna);
var7=netcdf.defVar(ncid,'ElevationReceiverAntenna','NC_FLOAT',dimid);
netcdf.putVar(ncid,var7,ElevationReceiverAntenna);
var8=netcdf.defVar(ncid,'ElevationTransmitterAntenna','NC_FLOAT',dimid);
netcdf.putVar(ncid,var8,ElevationTransmitterAntenna);
var9=netcdf.defVar(ncid,'SignalPower','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var9,SignalPower);
var10=netcdf.defVar(ncid,'NoisePower','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var10,NoisePower);
var11=netcdf.defVar(ncid,'BlackBodyPower','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var11,BlackBodyPower);
var12=netcdf.defVar(ncid,'SystemGainBB','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var12,SystemGainBB);
var13=netcdf.defVar(ncid,'SystemGainBBComp','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var13,SystemGainBBComp);
var14=netcdf.defVar(ncid,'RxTemperature','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var14,RxTemperature);
var15=netcdf.defVar(ncid,'InEclipse','NC_SHORT',dimid);
netcdf.putVar(ncid,var15,InEclipse);
var16=netcdf.defVar(ncid,'ReceiverPositionX','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var16,ReceiverPositionX);
var17=netcdf.defVar(ncid,'ReceiverPositionY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var17,ReceiverPositionY);
var18=netcdf.defVar(ncid,'ReceiverPositionZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var18,ReceiverPositionZ);
var19=netcdf.defVar(ncid,'ReceiverVelocityX','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var19,ReceiverVelocityX);
var20=netcdf.defVar(ncid,'ReceiverVelocityY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var20,ReceiverVelocityY);
var21=netcdf.defVar(ncid,'ReceiverVelocityZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var21,ReceiverVelocityZ);
var22=netcdf.defVar(ncid,'TransmitterPositionX','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var22,TransmitterPositionX);
var23=netcdf.defVar(ncid,'TransmitterPositionY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var23,TransmitterPositionY);
var24=netcdf.defVar(ncid,'TransmitterPositionZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var24,TransmitterPositionZ);
var25=netcdf.defVar(ncid,'TransmitterVelocityX','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var25,TransmitterVelocityX);
var26=netcdf.defVar(ncid,'TransmitterVelocityY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var26,TransmitterVelocityY);
var27=netcdf.defVar(ncid,'TransmitterVelocityZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var27,TransmitterVelocityZ);
var28=netcdf.defVar(ncid,'AttitudeRoll','NC_FLOAT',dimid);
netcdf.putVar(ncid,var28,AttitudeRoll);
var29=netcdf.defVar(ncid,'AttitudePitch','NC_FLOAT',dimid);
netcdf.putVar(ncid,var29,AttitudePitch);
var30=netcdf.defVar(ncid,'AttitudeYaw','NC_FLOAT',dimid);
netcdf.putVar(ncid,var30,AttitudeYaw);
var32=netcdf.defVar(ncid,'SPAzimuthARF','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var32,SPAzimuthARF);
var33=netcdf.defVar(ncid,'SPElevationARF','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var33,SPElevationARF);
var34=netcdf.defVar(ncid,'InternalTemperatures0','NC_FLOAT',dimid);
netcdf.putVar(ncid,var34,InternalTemperatures0);
var35=netcdf.defVar(ncid,'InternalTemperatures1','NC_FLOAT',dimid);
netcdf.putVar(ncid,var35,InternalTemperatures1);
var36=netcdf.defVar(ncid,'InternalTemperatures2','NC_FLOAT',dimid);
netcdf.putVar(ncid,var36,InternalTemperatures2);
var37=netcdf.defVar(ncid,'InternalTemperatures3','NC_FLOAT',dimid);
netcdf.putVar(ncid,var37,InternalTemperatures3);
var38=netcdf.defVar(ncid,'InternalTemperatures4','NC_FLOAT',dimid);
netcdf.putVar(ncid,var38,InternalTemperatures4);
var39=netcdf.defVar(ncid,'InternalTemperatures5','NC_FLOAT',dimid);
netcdf.putVar(ncid,var39,InternalTemperatures5);
var40=netcdf.defVar(ncid,'InternalTemperatures6','NC_FLOAT',dimid);
netcdf.putVar(ncid,var40,InternalTemperatures6);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'Name','HydroGNSS L1a directSignalPower File');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceConstellation','HYDROGNSS');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceNumber',SourceNumber);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'L1a_ProcessorVersion',L1a_ProcessorVersion);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'StartTime',string(StartTime));
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'StopTime',string(StopTime));


netcdf.close(ncid);

 end