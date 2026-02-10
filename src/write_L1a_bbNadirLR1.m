%% ************************************************************************
% MODULE NAME:      write_L1a_blackbodynadir.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      Aneesha Musunuri, N. Pierdicca
% @copyright:   Sapienza University of Rome
% Original version:      Sep 2021 - Sep 2022 by Aneesha Musunuri, N. Pierdicca
% Main updates:          Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Very Minor updates:    Jul 2023 by Davide Comite
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module writes the blackbodyNadir.nc L1a product of the Hydrognss
% E2ES.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function write_L1a_bbNadirLR1(IntegrationMidPointTime, geoSYSp, Rx, DDMOutputNumericalScaling_NLblackbody_DOWNmax, DDMPixelValueNoiseWholeDDM_down,...
                           DDM_NLBB_down, PRN, outpathname, foldername)

%% variable declaration
SourceNumber                         = uint32(1);
L1a_ProcessorVersion                 = string(1);
datecreated                          = datetime('now','Format','yyyy-MM-dd''T''hh:mm:ss');
ConfigurationNumber                  = ones(size(IntegrationMidPointTime));
ReflectionChannelNumber              = repmat(4,size(IntegrationMidPointTime));
FrontEnd                             = zeros(size(IntegrationMidPointTime));
SamplingFrequency                    = repmat(geoSYSp.SamplingFrequency_Hz,size(IntegrationMidPointTime));
FrontEndGainMode                     = zeros(size(IntegrationMidPointTime));
FrontEndGainControlSetting           = repmat(36,size(IntegrationMidPointTime));
CoherentIntegrationTime              = repmat(geoSYSp.cohIntTime,size(IntegrationMidPointTime));
CodeDelaySpacingSamplesBetweenPixels = repmat(uint16(geoSYSp.CodeDelaySpacingSamplesBetweenPixels(1)),size(IntegrationMidPointTime));
IncoherentIntegrations               = repmat((geoSYSp.incoIntTime),size(IntegrationMidPointTime));
NumberOfDelayPixels                  = int32(repmat(geoSYSp.nx_Rdelay,size(IntegrationMidPointTime)));
NumberOfDopplerPixels                = int32(repmat(geoSYSp.ny_Dfreq,size(IntegrationMidPointTime)));
DopplerResolution                    = repmat(geoSYSp.fpas,size(IntegrationMidPointTime));
% DelayResolution=repmat((10^3*geoSYSp.Dly_RNGmax/(geoSYSp.nx_Rdelay)),size(IntegrationMidPointTime));
% DelayResolution=repmat((double(double(CodeDelaySpacingSamplesBetweenPixels)/SamplingFrequency*1e9)),size(IntegrationMidPointTime));
DelayResolution                      = repmat(double(10^3*geoSYSp.mpas(1)),size(IntegrationMidPointTime));
TrackingOffsetDelayInPixels          = uint32(repmat(geoSYSp.delayOffset_bin,size(IntegrationMidPointTime)));
TrackingOffsetDopplerHz              = zeros(size(IntegrationMidPointTime));
TrackingOffsetDopplerInPixels        = zeros(size(IntegrationMidPointTime));
TrackingOffsetDelayNs                = repmat(10^3*geoSYSp.delayOffset_microsec(1),size(IntegrationMidPointTime));
[l,~]                                = size(IntegrationMidPointTime);
ADCPercentagePositive                = zeros(size(IntegrationMidPointTime));
ADCOffset                            = repmat(0.039,size(IntegrationMidPointTime));
ADCPercentageMagnitude               = repmat(29.4961,size(IntegrationMidPointTime));
ReceiverPositionX   = Rx(:,1);
ReceiverPositionY   = Rx(:,2);
ReceiverPositionZ   = Rx(:,3);
InternalTemperature = geoSYSp.TRx_K - 273.15;
[InternalTemperatures0, InternalTemperatures1, InternalTemperatures2, InternalTemperatures3, InternalTemperatures4, InternalTemperatures5, InternalTemperatures6]=...
 deal(repmat(InternalTemperature, (size(ReceiverPositionZ))));


DDMOutputNumericalScaling  = DDMOutputNumericalScaling_NLblackbody_DOWNmax;
DDMPixelValueNoiseWholeDDM = DDMPixelValueNoiseWholeDDM_down;
maxddm                     = max(max(DDM_NLBB_down));
ddm1                       = zeros(size(DDM_NLBB_down));

for ki=1:l
    ddm1(:,:,ki) = ((DDM_NLBB_down(:,:,ki).*65535./(maxddm(ki))));
    % ddm(:,:,ki)=DDM_a(:,:,ki) .*65535 ./ (max(max(DDM_a(:,:,ki))));
end

DDM     = uint16(ddm1);
Delay   = (0:geoSYSp.nx_Rdelay-1)';
[a1,~]  = size(Delay);
Doppler = (-fix(geoSYSp.ny_Dfreq/2) : round(geoSYSp.ny_Dfreq/2)-1)';
[a2,~]  = size(Doppler);
signalpolarization = 'LHCP';
filename = ['blackbodyNadirL1E1' signalpolarization '.nc'];

% filename = 
%% creating the file
ncid  = netcdf.create(strjoin([outpathname '\' filename],''),'netcdf4');
[a,~] = size(IntegrationMidPointTime);
dimid = netcdf.defDim(ncid,'DateTime',a);

%% variables
var1=netcdf.defVar(ncid,'IntegrationMidPointTime','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var1,IntegrationMidPointTime);
var2=netcdf.defVar(ncid,'ConfigurationNumber','NC_INT',dimid);
netcdf.putVar(ncid,var2,ConfigurationNumber);
var3=netcdf.defVar(ncid,'ReflectionChannelNumber','NC_SHORT',dimid);
netcdf.putVar(ncid,var3,ReflectionChannelNumber);
var4=netcdf.defVar(ncid,'FrontEnd','NC_SHORT',dimid);
netcdf.putVar(ncid,var4,FrontEnd);
var5=netcdf.defVar(ncid,'SamplingFrequency','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var5,SamplingFrequency);
var6=netcdf.defVar(ncid,'FrontEndGainMode','NC_SHORT',dimid);
netcdf.putVar(ncid,var6,FrontEndGainMode);
var7=netcdf.defVar(ncid,'FrontEndGainControlSetting','NC_SHORT',dimid);
netcdf.putVar(ncid,var7,FrontEndGainControlSetting);
var8=netcdf.defVar(ncid,'CoherentIntegrationTime','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var8,CoherentIntegrationTime);
var9=netcdf.defVar(ncid,'DelayResolution','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var9,DelayResolution);
var10=netcdf.defVar(ncid,'CodeDelaySpacingSamplesBetweenPixels','NC_SHORT',dimid);
netcdf.putVar(ncid,var10,CodeDelaySpacingSamplesBetweenPixels);
var11=netcdf.defVar(ncid,'DopplerResolution','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var11,DopplerResolution);
var12=netcdf.defVar(ncid,'IncoherentIntegrations','NC_INT',dimid);
netcdf.putVar(ncid,var12,IncoherentIntegrations);
var13=netcdf.defVar(ncid,'NumberOfDelayPixels','NC_INT',dimid);
netcdf.putVar(ncid,var13,NumberOfDelayPixels);
var14=netcdf.defVar(ncid,'NumberOfDopplerPixels','NC_INT',dimid);
netcdf.putVar(ncid,var14,NumberOfDopplerPixels);
var15=netcdf.defVar(ncid,'TrackingOffsetDelayNs','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var15,TrackingOffsetDelayNs);
var16=netcdf.defVar(ncid,'TrackingOffsetDelayInPixels','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var16,TrackingOffsetDelayInPixels);
var17=netcdf.defVar(ncid,'TrackingOffsetDopplerHz','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var17,TrackingOffsetDopplerHz);
var18=netcdf.defVar(ncid,'TrackingOffsetDopplerInPixels','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var18,TrackingOffsetDopplerInPixels);
var19=netcdf.defVar(ncid,'PRN','NC_USHORT',dimid);
netcdf.putVar(ncid,var19,PRN);

var21=netcdf.defVar(ncid,'ADCPercentagePositive','NC_FLOAT',dimid);
netcdf.putVar(ncid,var21,ADCPercentagePositive);
var22=netcdf.defVar(ncid,'ADCOffset','NC_FLOAT',dimid);
netcdf.putVar(ncid,var22,ADCOffset);
var23=netcdf.defVar(ncid,'ADCPercentageMagnitude','NC_FLOAT',dimid);
netcdf.putVar(ncid,var23,ADCPercentageMagnitude);
var24=netcdf.defVar(ncid,'DDMOutputNumericalScaling','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var24,DDMOutputNumericalScaling);
var25=netcdf.defVar(ncid,'DDMPixelValueNoiseWholeDDM','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var25,DDMPixelValueNoiseWholeDDM);
var26=netcdf.defVar(ncid,'InternalTemperatures0','NC_FLOAT',dimid);
netcdf.putVar(ncid,var26,InternalTemperatures0);
var27=netcdf.defVar(ncid,'InternalTemperatures1','NC_FLOAT',dimid);
netcdf.putVar(ncid,var27,InternalTemperatures1);
var28=netcdf.defVar(ncid,'InternalTemperatures2','NC_FLOAT',dimid);
netcdf.putVar(ncid,var28,InternalTemperatures2);
var29=netcdf.defVar(ncid,'InternalTemperatures3','NC_FLOAT',dimid);
netcdf.putVar(ncid,var29,InternalTemperatures3);
var29_a=netcdf.defVar(ncid,'InternalTemperatures4','NC_FLOAT',dimid);
netcdf.putVar(ncid,var29_a,InternalTemperatures4);
var29_b=netcdf.defVar(ncid,'InternalTemperatures5','NC_FLOAT',dimid);
netcdf.putVar(ncid,var29_b,InternalTemperatures5);
var29_c=netcdf.defVar(ncid,'InternalTemperatures6','NC_FLOAT',dimid);
netcdf.putVar(ncid,var29_c,InternalTemperatures6);
var30=netcdf.defVar(ncid,'ReceiverPositionX','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var30,ReceiverPositionX);
var31=netcdf.defVar(ncid,'ReceiverPositionY','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var31,ReceiverPositionY);
var32=netcdf.defVar(ncid,'ReceiverPositionZ','NC_DOUBLE',dimid);
netcdf.putVar(ncid,var32,ReceiverPositionZ);
dimid2 = netcdf.defDim(ncid,'Delay',a1);
var33=netcdf.defVar(ncid,'Delay','NC_DOUBLE',dimid2);
netcdf.putVar(ncid,var33,Delay);
dimid3 = netcdf.defDim(ncid,'Doppler',a2);
var34=netcdf.defVar(ncid,'Doppler','NC_DOUBLE',dimid3);
netcdf.putVar(ncid,var34,Doppler);
dimid4=[dimid2 dimid3 dimid];
var35=netcdf.defVar(ncid,'DDM','NC_USHORT',dimid4);
netcdf.putVar(ncid,var35,DDM);

%% attributes
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'Name','HydroGNSS L1a blackbodyNadirL1E1LHCP.nc');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'date_created',string(datecreated));
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'License','HydroGNSS GNSS-R Dataset');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'FileIDCode',foldername);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceConstellation','HYDROGNSS');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'SourceNumber',SourceNumber);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'L1a_ProcessorVersion',L1a_ProcessorVersion);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'BlackbodyType',uint32(1));
netcdf.close(ncid);

 end