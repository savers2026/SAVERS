%% ************************************************************************
% MODULE NAME:      readConfigurationFile.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome
% Original version:      Dec 2022 by Laura Dente
% Main updates:          Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module reads the HSAVERS configuration files and return the
% paraemters to the main routine.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [configParam,instrumentConfParam] = readConfigurationFile(HSAVERSConfigurationFilename,InstrumentConfigurationFilename)

%% read HSAVERSConfiguration file
confFile_table=readtable(HSAVERSConfigurationFilename,'Format','auto','ReadRowNames',true);
confFile_struct=table2struct(confFile_table,'ToScalar',true);
confFile_struct.RowNames = confFile_table.Properties.RowNames;

indexBioGeoFact = find(strcmpi(confFile_struct.RowNames,{'biogeofact'}));
configParam.biogeofact=str2double(confFile_struct.VALUE(indexBioGeoFact));

indexChadDEM3 = find(strcmpi(confFile_struct.RowNames,{'Chad.DEM3arcsec'}));
configParam.Chad.DEM3arcsec=char(confFile_struct.VALUE(indexChadDEM3));
indexChadDEM9 = find(strcmpi(confFile_struct.RowNames,{'Chad.DEM9arcsec'}));
configParam.Chad.DEM9arcsec=char(confFile_struct.VALUE(indexChadDEM9));
indexChadDEM30 = find(strcmpi(confFile_struct.RowNames,{'Chad.DEM30arcsec'}));
configParam.Chad.DEM30arcsec=char(confFile_struct.VALUE(indexChadDEM30));

indexChadMinLat = find(strcmpi(confFile_struct.RowNames,{'Chad.SPcoordLimits.minLat'}));
ChadminLat=str2double(confFile_struct.VALUE(indexChadMinLat));
indexChadMaxLat = find(strcmpi(confFile_struct.RowNames,{'Chad.SPcoordLimits.maxLat'}));
ChadmaxLat=str2double(confFile_struct.VALUE(indexChadMaxLat));
indexChadMinLon = find(strcmpi(confFile_struct.RowNames,{'Chad.SPcoordLimits.minLon'}));
ChadminLon=str2double(confFile_struct.VALUE(indexChadMinLon));
indexChadMaxLon = find(strcmpi(confFile_struct.RowNames,{'Chad.SPcoordLimits.maxLon'}));
ChadmaxLon=str2double(confFile_struct.VALUE(indexChadMaxLon));
configParam.Chad.SPcoordLimits=[ChadminLat ChadmaxLat ChadminLon ChadmaxLon];

indexCherskyDEM3 = find(strcmpi(confFile_struct.RowNames,{'Chersky.DEM3arcsec'}));
configParam.Chersky.DEM3arcsec=char(confFile_struct.VALUE(indexCherskyDEM3));
indexCherskyDEM9 = find(strcmpi(confFile_struct.RowNames,{'Chersky.DEM9arcsec'}));
configParam.Chersky.DEM9arcsec=char(confFile_struct.VALUE(indexCherskyDEM9));
indexCherskyDEM30 = find(strcmpi(confFile_struct.RowNames,{'Chersky.DEM30arcsec'}));
configParam.Chersky.DEM30arcsec=char(confFile_struct.VALUE(indexCherskyDEM30));

indexCherskyMinLat = find(strcmpi(confFile_struct.RowNames,{'Chersky.SPcoordLimits.minLat'}));
CherskyminLat=str2double(confFile_struct.VALUE(indexCherskyMinLat));
indexCherskyMaxLat = find(strcmpi(confFile_struct.RowNames,{'Chersky.SPcoordLimits.maxLat'}));
CherskymaxLat=str2double(confFile_struct.VALUE(indexCherskyMaxLat));
indexCherskyMinLon = find(strcmpi(confFile_struct.RowNames,{'Chersky.SPcoordLimits.minLon'}));
CherskyminLon=str2double(confFile_struct.VALUE(indexCherskyMinLon));
indexCherskyMaxLon = find(strcmpi(confFile_struct.RowNames,{'Chersky.SPcoordLimits.maxLon'}));
CherskymaxLon=str2double(confFile_struct.VALUE(indexCherskyMaxLon));
configParam.Chersky.SPcoordLimits=[CherskyminLat CherskymaxLat CherskyminLon CherskymaxLon];

indexCodeRepetitionFreq = find(strcmpi(confFile_struct.RowNames,{'CodeRepetitionFreq_Hz'}));
configParam.CodeRepetitionFreq_Hz=str2double(confFile_struct.VALUE(indexCodeRepetitionFreq));

indexdang = find(strcmpi(confFile_struct.RowNames,{'dang'}));
configParam.dang=str2double(confFile_struct.VALUE(indexdang));

indexdelayOffsetFactor = find(strcmpi(confFile_struct.RowNames,{'delayOffsetFactor'}));
configParam.delayOffsetFactor=str2double(confFile_struct.VALUE(indexdelayOffsetFactor));

indexEvergladesDEM3 = find(strcmpi(confFile_struct.RowNames,{'Everglades.DEM3arcsec'}));
configParam.Everglades.DEM3arcsec=char(confFile_struct.VALUE(indexEvergladesDEM3));
indexEvergladesDEM9 = find(strcmpi(confFile_struct.RowNames,{'Everglades.DEM9arcsec'}));
configParam.Everglades.DEM9arcsec=char(confFile_struct.VALUE(indexEvergladesDEM9));
indexEvergladesDEM30 = find(strcmpi(confFile_struct.RowNames,{'Everglades.DEM30arcsec'}));
configParam.Everglades.DEM30arcsec=char(confFile_struct.VALUE(indexEvergladesDEM30));

indexEvergladesMinLat = find(strcmpi(confFile_struct.RowNames,{'Everglades.SPcoordLimits.minLat'}));
EvergladesminLat=str2double(confFile_struct.VALUE(indexEvergladesMinLat));
indexEvergladesMaxLat = find(strcmpi(confFile_struct.RowNames,{'Everglades.SPcoordLimits.maxLat'}));
EvergladesmaxLat=str2double(confFile_struct.VALUE(indexEvergladesMaxLat));
indexEvergladesMinLon = find(strcmpi(confFile_struct.RowNames,{'Everglades.SPcoordLimits.minLon'}));
EvergladesminLon=str2double(confFile_struct.VALUE(indexEvergladesMinLon));
indexEvergladesMaxLon = find(strcmpi(confFile_struct.RowNames,{'Everglades.SPcoordLimits.maxLon'}));
EvergladesmaxLon=str2double(confFile_struct.VALUE(indexEvergladesMaxLon));
configParam.Everglades.SPcoordLimits=[EvergladesminLat EvergladesmaxLat EvergladesminLon EvergladesmaxLon];

indexfilenameCoverMap = find(strcmpi(confFile_struct.RowNames,{'filenameCoverMap'}));
configParam.filenameCoverMap=char(confFile_struct.VALUE(indexfilenameCoverMap));

indexfilenameDEM5km = find(strcmpi(confFile_struct.RowNames,{'filenameDEM5km'}));
configParam.filenameDEM5km=char(confFile_struct.VALUE(indexfilenameDEM5km));

indexfilenameIonisationLevelParam = find(strcmpi(confFile_struct.RowNames,{'filenameIonisationLevelParam'}));
configParam.filenameIonisationLevelParam=char(confFile_struct.VALUE(indexfilenameIonisationLevelParam));

indexFinlandDEM3 = find(strcmpi(confFile_struct.RowNames,{'Finland.DEM3arcsec'}));
configParam.Finland.DEM3arcsec=char(confFile_struct.VALUE(indexFinlandDEM3));
indexFinlandDEM9 = find(strcmpi(confFile_struct.RowNames,{'Finland.DEM9arcsec'}));
configParam.Finland.DEM9arcsec=char(confFile_struct.VALUE(indexFinlandDEM9));
indexFinlandDEM30 = find(strcmpi(confFile_struct.RowNames,{'Finland.DEM30arcsec'}));
configParam.Finland.DEM30arcsec=char(confFile_struct.VALUE(indexFinlandDEM30));

indexFinlandMinLat = find(strcmpi(confFile_struct.RowNames,{'Finland.SPcoordLimits.minLat'}));
FinlandminLat=str2double(confFile_struct.VALUE(indexFinlandMinLat));
indexFinlandMaxLat = find(strcmpi(confFile_struct.RowNames,{'Finland.SPcoordLimits.maxLat'}));
FinlandmaxLat=str2double(confFile_struct.VALUE(indexFinlandMaxLat));
indexFinlandMinLon = find(strcmpi(confFile_struct.RowNames,{'Finland.SPcoordLimits.minLon'}));
FinlandminLon=str2double(confFile_struct.VALUE(indexFinlandMinLon));
indexFinlandMaxLon = find(strcmpi(confFile_struct.RowNames,{'Finland.SPcoordLimits.maxLon'}));
FinlandmaxLon=str2double(confFile_struct.VALUE(indexFinlandMaxLon));
configParam.Finland.SPcoordLimits=[FinlandminLat FinlandmaxLat FinlandminLon FinlandmaxLon];

indexfrequencyGPSL1_GalileoE1c = find(strcmpi(confFile_struct.RowNames,{'frequencyGPSL1_GalileoE1c_GHz'}));
configParam.frequencyGPSL1_GalileoE1c_GHz=str2double(confFile_struct.VALUE(indexfrequencyGPSL1_GalileoE1c));

indexfrequencyGPSL5_GalileoE5a = find(strcmpi(confFile_struct.RowNames,{'frequencyGPSL5_GalileoE5a_GHz'}));
configParam.frequencyGPSL5_GalileoE5a_GHz=str2double(confFile_struct.VALUE(indexfrequencyGPSL5_GalileoE5a));

indexGabonDEM3 = find(strcmpi(confFile_struct.RowNames,{'Gabon.DEM3arcsec'}));
configParam.Gabon.DEM3arcsec=char(confFile_struct.VALUE(indexGabonDEM3));
indexGabonDEM9 = find(strcmpi(confFile_struct.RowNames,{'Gabon.DEM9arcsec'}));
configParam.Gabon.DEM9arcsec=char(confFile_struct.VALUE(indexGabonDEM9));
indexGabonDEM30 = find(strcmpi(confFile_struct.RowNames,{'Gabon.DEM30arcsec'}));
configParam.Gabon.DEM30arcsec=char(confFile_struct.VALUE(indexGabonDEM30));

indexGabonMinLat = find(strcmpi(confFile_struct.RowNames,{'Gabon.SPcoordLimits.minLat'}));
GabonminLat=str2double(confFile_struct.VALUE(indexGabonMinLat));
indexGabonMaxLat = find(strcmpi(confFile_struct.RowNames,{'Gabon.SPcoordLimits.maxLat'}));
GabonmaxLat=str2double(confFile_struct.VALUE(indexGabonMaxLat));
indexGabonMinLon = find(strcmpi(confFile_struct.RowNames,{'Gabon.SPcoordLimits.minLon'}));
GabonminLon=str2double(confFile_struct.VALUE(indexGabonMinLon));
indexGabonMaxLon = find(strcmpi(confFile_struct.RowNames,{'Gabon.SPcoordLimits.maxLon'}));
GabonmaxLon=str2double(confFile_struct.VALUE(indexGabonMaxLon));
configParam.Gabon.SPcoordLimits=[GabonminLat GabonmaxLat GabonminLon GabonmaxLon];

indexgeoboundedWindowWidthInDelay = find(strcmpi(confFile_struct.RowNames,{'geoboundedWindowWidthInDelay'}));
geoboundedWindowWidthInDelay=str2double(confFile_struct.VALUE(indexgeoboundedWindowWidthInDelay));
indexgeoboundedWindowWidthInFreq = find(strcmpi(confFile_struct.RowNames,{'geoboundedWindowWidthInFreq'}));
geoboundedWindowWidthInFreq=str2double(confFile_struct.VALUE(indexgeoboundedWindowWidthInFreq));
configParam.geoboundedWindowWidth = [geoboundedWindowWidthInDelay geoboundedWindowWidthInFreq];

indexRxEuler_phi_angle = find(strcmpi(confFile_struct.RowNames,{'RxEuler_phi_angle'}));
RxEuler_phi_angle=str2double(confFile_struct.VALUE(indexRxEuler_phi_angle));
indexRxEuler_theta_angle = find(strcmpi(confFile_struct.RowNames,{'RxEuler_theta_angle'}));
RxEuler_theta_angle=str2double(confFile_struct.VALUE(indexRxEuler_theta_angle));
indexRxEuler_psi_angle = find(strcmpi(confFile_struct.RowNames,{'RxEuler_psi_angle'}));
RxEuler_psi_angle=str2double(confFile_struct.VALUE(indexRxEuler_psi_angle));
configParam.RxEuler_angles=[RxEuler_phi_angle RxEuler_theta_angle RxEuler_psi_angle];

indexSamplingFrequency = find(strcmpi(confFile_struct.RowNames,{'SamplingFrequency_Hz'}));
configParam.SamplingFrequency_Hz=str2double(confFile_struct.VALUE(indexSamplingFrequency));

indexsizeInputMapWindow = find(strcmpi(confFile_struct.RowNames,{'sizeInputMapWindow'}));
configParam.sizeInputMapWindow=str2double(confFile_struct.VALUE(indexsizeInputMapWindow));

indexsizeInputMapWindowCoarseMode = find(strcmpi(confFile_struct.RowNames,{'sizeInputMapWindowCoarseMode'}));
configParam.sizeInputMapWindowCoarseMode=str2double(confFile_struct.VALUE(indexsizeInputMapWindowCoarseMode));

indexsizeInputMapWindowHighLatitudes = find(strcmpi(confFile_struct.RowNames,{'sizeInputMapWindowHighLatitudes'}));
configParam.sizeInputMapWindowHighLatitudes=str2double(confFile_struct.VALUE(indexsizeInputMapWindowHighLatitudes));

indexsizeInputMapWindowHighLatitudesCoarseMode = find(strcmpi(confFile_struct.RowNames,{'sizeInputMapWindowHighLatitudesCoarseMode'}));
configParam.sizeInputMapWindowHighLatitudesCoarseMode=str2double(confFile_struct.VALUE(indexsizeInputMapWindowHighLatitudesCoarseMode));

indexsizeIntegrationAreaCoarseMode = find(strcmpi(confFile_struct.RowNames,{'sizeIntegrationAreaCoarseMode'}));
configParam.sizeIntegrationAreaCoarseMode=str2double(confFile_struct.VALUE(indexsizeIntegrationAreaCoarseMode));

indexSLVDEM3 = find(strcmpi(confFile_struct.RowNames,{'SLV.DEM3arcsec'}));
configParam.SLV.DEM3arcsec=char(confFile_struct.VALUE(indexSLVDEM3));
indexSLVDEM9 = find(strcmpi(confFile_struct.RowNames,{'SLV.DEM9arcsec'}));
configParam.SLV.DEM9arcsec=char(confFile_struct.VALUE(indexSLVDEM9));
indexSLVDEM30 = find(strcmpi(confFile_struct.RowNames,{'SLV.DEM30arcsec'}));
configParam.SLV.DEM30arcsec=char(confFile_struct.VALUE(indexSLVDEM30));

indexSLVMinLat = find(strcmpi(confFile_struct.RowNames,{'SLV.SPcoordLimits.minLat'}));
SLVminLat=str2double(confFile_struct.VALUE(indexSLVMinLat));
indexSLVMaxLat = find(strcmpi(confFile_struct.RowNames,{'SLV.SPcoordLimits.maxLat'}));
SLVmaxLat=str2double(confFile_struct.VALUE(indexSLVMaxLat));
indexSLVMinLon = find(strcmpi(confFile_struct.RowNames,{'SLV.SPcoordLimits.minLon'}));
SLVminLon=str2double(confFile_struct.VALUE(indexSLVMinLon));
indexSLVMaxLon = find(strcmpi(confFile_struct.RowNames,{'SLV.SPcoordLimits.maxLon'}));
SLVmaxLon=str2double(confFile_struct.VALUE(indexSLVMaxLon));
configParam.SLV.SPcoordLimits=[SLVminLat SLVmaxLat SLVminLon SLVmaxLon];

indexSolarActivityLevel = find(strcmpi(confFile_struct.RowNames,{'SolarActivityLevel'}));
configParam.SolarActivityLevel=char(confFile_struct.VALUE(indexSolarActivityLevel));

indexTx_nomRollAngle = find(strcmpi(confFile_struct.RowNames,{'Tx_nomRollAngle'}));
Tx_nomRollAngle=str2double(confFile_struct.VALUE(indexTx_nomRollAngle));
indexTx_nomPitchAngle = find(strcmpi(confFile_struct.RowNames,{'Tx_nomPitchAngle'}));
Tx_nomPitchAngle=str2double(confFile_struct.VALUE(indexTx_nomPitchAngle));
indexTx_nomYawAngle = find(strcmpi(confFile_struct.RowNames,{'Tx_nomYawAngle'}));
Tx_nomYawAngle=str2double(confFile_struct.VALUE(indexTx_nomYawAngle));
configParam.Tx_nomAtt_angles=[Tx_nomRollAngle Tx_nomPitchAngle Tx_nomYawAngle];

indexTxAntennaPattern_E1E5 = find(strcmpi(confFile_struct.RowNames,{'TxAntennaPattern_E1_E5'}));
configParam.TxAntennaPattern_E1E5=char(confFile_struct.VALUE(indexTxAntennaPattern_E1E5));

indexTxAntennaPattern_L1L5 = find(strcmpi(confFile_struct.RowNames,{'TxAntennaPattern_L1_L5'}));
configParam.TxAntennaPattern_L1L5=char(confFile_struct.VALUE(indexTxAntennaPattern_L1L5));

indexTxPower_E1bc_dBW = find(strcmpi(confFile_struct.RowNames,{'TxPower_E1bc_dBW'}));
configParam.TxPower_E1bc_dBW=str2double(confFile_struct.VALUE(indexTxPower_E1bc_dBW));

indexTxPower_E5a_dBW = find(strcmpi(confFile_struct.RowNames,{'TxPower_E5a_dBW'}));
configParam.TxPower_E5a_dBW=str2double(confFile_struct.VALUE(indexTxPower_E5a_dBW));

indexTxPower_L1_dBW = find(strcmpi(confFile_struct.RowNames,{'TxPower_L1_dBW'}));
configParam.TxPower_L1_dBW=str2double(confFile_struct.VALUE(indexTxPower_L1_dBW));

indexTxPower_L5_dBW_temp = find(strcmpi(confFile_struct.RowNames,{'TxPower_L5_dBW_temp'}));
configParam.TxPower_L5_dBW_temp=str2double(confFile_struct.VALUE(indexTxPower_L5_dBW_temp));

indexYucatanDEM3 = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM3arcsec'}));
configParam.Yucatan.DEM3arcsec=char(confFile_struct.VALUE(indexYucatanDEM3));
indexYucatanDEM9 = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM9arcsec'}));
configParam.Yucatan.DEM9arcsec=char(confFile_struct.VALUE(indexYucatanDEM9));
indexYucatanDEM30 = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM30arcsec'}));
configParam.Yucatan.DEM30arcsec=char(confFile_struct.VALUE(indexYucatanDEM30));
indexYucatanDEM3flooded = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM3arcsec_flooded'}));
configParam.Yucatan.DEM3arcsec_flooded=char(confFile_struct.VALUE(indexYucatanDEM3flooded));
indexYucatanDEM9flooded = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM9arcsec_flooded'}));
configParam.Yucatan.DEM9arcsec_flooded=char(confFile_struct.VALUE(indexYucatanDEM9flooded));
indexYucatanDEM30flooded = find(strcmpi(confFile_struct.RowNames,{'Yucatan.DEM30arcsec_flooded'}));
configParam.Yucatan.DEM30arcsec_flooded=char(confFile_struct.VALUE(indexYucatanDEM30flooded));

indexYucatanMinLat = find(strcmpi(confFile_struct.RowNames,{'Yucatan.SPcoordLimits.minLat'}));
YucatanminLat=str2double(confFile_struct.VALUE(indexYucatanMinLat));
indexYucatanMaxLat = find(strcmpi(confFile_struct.RowNames,{'Yucatan.SPcoordLimits.maxLat'}));
YucatanmaxLat=str2double(confFile_struct.VALUE(indexYucatanMaxLat));
indexYucatanMinLon = find(strcmpi(confFile_struct.RowNames,{'Yucatan.SPcoordLimits.minLon'}));
YucatanminLon=str2double(confFile_struct.VALUE(indexYucatanMinLon));
indexYucatanMaxLon = find(strcmpi(confFile_struct.RowNames,{'Yucatan.SPcoordLimits.maxLon'}));
YucatanmaxLon=str2double(confFile_struct.VALUE(indexYucatanMaxLon));
configParam.Yucatan.SPcoordLimits=[YucatanminLat YucatanmaxLat YucatanminLon YucatanmaxLon];

%% read InstrumentConfiguration file
optsReadTable_instFile=detectImportOptions(InstrumentConfigurationFilename);
% L1/E1
optsReadTable_instFile.Sheet = 'parameters_L1E1';
instFile_table=readtable(InstrumentConfigurationFilename, optsReadTable_instFile,'ReadRowNames',true,"UseExcel", false);
instFile_struct=table2struct(instFile_table,'ToScalar',true);
instFile_struct.RowNames = instFile_table.Properties.RowNames;

indexF_dB_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'F_dB_LHCP_1'}));
F_dB_LHCP_1 = instFile_struct.VALUE(indexF_dB_LHCP_1);
indexF_dB_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'F_dB_RHCP_1'}));
F_dB_RHCP_1 = instFile_struct.VALUE(indexF_dB_RHCP_1);

indexgG_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'gG_LHCP_1'}));
gG_LHCP_1=instFile_struct.VALUE(indexgG_LHCP_1);
indexgG_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'gG_RHCP_1'}));
gG_RHCP_1=instFile_struct.VALUE(indexgG_RHCP_1);
indexgG_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'gG_RHCP_1_zenith'}));
gG_RHCP_1_zenith=instFile_struct.VALUE(indexgG_RHCP_1_zenith);

indexgNF_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'gNF_LHCP_1'}));
gNF_LHCP_1 = instFile_struct.VALUE(indexgNF_LHCP_1);
indexgNF_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'gNF_RHCP_1'}));
gNF_RHCP_1 = instFile_struct.VALUE(indexgNF_RHCP_1);
indexgNF_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'gNF_RHCP_1_zenith'}));
gNF_RHCP_1_zenith = instFile_struct.VALUE(indexgNF_RHCP_1_zenith);

indexGsys_dB_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_LHCP_1'}));
Gsys_dB_LHCP_1 = instFile_struct.VALUE(indexGsys_dB_LHCP_1);
indexGsys_dB_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_RHCP_1'}));
Gsys_dB_RHCP_1 = instFile_struct.VALUE(indexGsys_dB_RHCP_1);
indexGsys_dB_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_RHCP_1_zenith'}));
Gsys_dB_RHCP_1_zenith = instFile_struct.VALUE(indexGsys_dB_RHCP_1_zenith);

indexLN_dB_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'LN_dB_LHCP_1'}));
LN_dB_LHCP_1=instFile_struct.VALUE(indexLN_dB_LHCP_1);
indexLN_dB_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'LN_dB_RHCP_1'}));
LN_dB_RHCP_1=instFile_struct.VALUE(indexLN_dB_RHCP_1);
indexLN_dB_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'LN_dB_RHCP_1_zenith'}));
LN_dB_RHCP_1_zenith=instFile_struct.VALUE(indexLN_dB_RHCP_1_zenith);

indexLS_dB_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'LS_dB_LHCP_1'}));
LS_dB_LHCP_1=instFile_struct.VALUE(indexLS_dB_LHCP_1);
indexLS_dB_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'LS_dB_RHCP_1'}));
LS_dB_RHCP_1=instFile_struct.VALUE(indexLS_dB_RHCP_1);
indexLS_dB_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'LS_dB_RHCP_1_zenith'}));
LS_dB_RHCP_1_zenith=instFile_struct.VALUE(indexLS_dB_RHCP_1_zenith);

indexmG_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'mG_LHCP_1'}));
mG_LHCP_1 = instFile_struct.VALUE(indexmG_LHCP_1);
indexmG_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'mG_RHCP_1'}));
mG_RHCP_1 = instFile_struct.VALUE(indexmG_RHCP_1);
indexmG_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'mG_RHCP_1_zenith'}));
mG_RHCP_1_zenith = instFile_struct.VALUE(indexmG_RHCP_1_zenith);

indexmNF_LHCP_1 = find(strcmpi(instFile_struct.RowNames,{'mNF_LHCP_1'}));
mNF_LHCP_1=instFile_struct.VALUE(indexmNF_LHCP_1);
indexmNF_RHCP_1 = find(strcmpi(instFile_struct.RowNames,{'mNF_RHCP_1'}));
mNF_RHCP_1=instFile_struct.VALUE(indexmNF_RHCP_1);
indexmNF_RHCP_1_zenith = find(strcmpi(instFile_struct.RowNames,{'mNF_RHCP_1_zenith'}));
mNF_RHCP_1_zenith=instFile_struct.VALUE(indexmNF_RHCP_1_zenith);

indexTL_K = find(strcmpi(instFile_struct.RowNames,{'TL_K'}));
instrumentConfParam.TL_K=instFile_struct.VALUE(indexTL_K);

indexTref_K = find(strcmpi(instFile_struct.RowNames,{'Tref_K'}));
instrumentConfParam.Tref_K=instFile_struct.VALUE(indexTref_K);

indexTRx_K = find(strcmpi(instFile_struct.RowNames,{'TRx_K'}));
instrumentConfParam.TRx_K=instFile_struct.VALUE(indexTRx_K);

% L5/E5
optsReadTable_instFile.Sheet = 'parameters_L5E5';
instFile_table=readtable(InstrumentConfigurationFilename, optsReadTable_instFile,'ReadRowNames',true,"UseExcel", false);
instFile_struct=table2struct(instFile_table,'ToScalar',true);
instFile_struct.RowNames = instFile_table.Properties.RowNames;

indexF_dB_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'F_dB_LHCP_5'}));
F_dB_LHCP_5 = instFile_struct.VALUE(indexF_dB_LHCP_5);
indexF_dB_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'F_dB_RHCP_5'}));
F_dB_RHCP_5 = instFile_struct.VALUE(indexF_dB_RHCP_5);

indexgG_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'gG_LHCP_5'}));
gG_LHCP_5=instFile_struct.VALUE(indexgG_LHCP_5);
indexgG_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'gG_RHCP_5'}));
gG_RHCP_5=instFile_struct.VALUE(indexgG_RHCP_5);
indexgG_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'gG_RHCP_5_zenith'}));
gG_RHCP_5_zenith=instFile_struct.VALUE(indexgG_RHCP_5_zenith);

indexgNF_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'gNF_LHCP_5'}));
gNF_LHCP_5 = instFile_struct.VALUE(indexgNF_LHCP_5);
indexgNF_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'gNF_RHCP_5'}));
gNF_RHCP_5 = instFile_struct.VALUE(indexgNF_RHCP_5);
indexgNF_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'gNF_RHCP_5_zenith'}));
gNF_RHCP_5_zenith = instFile_struct.VALUE(indexgNF_RHCP_5_zenith);

indexGsys_dB_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_LHCP_5'}));
Gsys_dB_LHCP_5 = instFile_struct.VALUE(indexGsys_dB_LHCP_5);
indexGsys_dB_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_RHCP_5'}));
Gsys_dB_RHCP_5 = instFile_struct.VALUE(indexGsys_dB_RHCP_5);
indexGsys_dB_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'Gsys_dB_RHCP_5_zenith'}));
Gsys_dB_RHCP_5_zenith = instFile_struct.VALUE(indexGsys_dB_RHCP_5_zenith);

indexLN_dB_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'LN_dB_LHCP_5'}));
LN_dB_LHCP_5=instFile_struct.VALUE(indexLN_dB_LHCP_5);
indexLN_dB_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'LN_dB_RHCP_5'}));
LN_dB_RHCP_5=instFile_struct.VALUE(indexLN_dB_RHCP_5);
indexLN_dB_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'LN_dB_RHCP_5_zenith'}));
LN_dB_RHCP_5_zenith=instFile_struct.VALUE(indexLN_dB_RHCP_5_zenith);

indexLS_dB_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'LS_dB_LHCP_5'}));
LS_dB_LHCP_5=instFile_struct.VALUE(indexLS_dB_LHCP_5);
indexLS_dB_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'LS_dB_RHCP_5'}));
LS_dB_RHCP_5=instFile_struct.VALUE(indexLS_dB_RHCP_5);
indexLS_dB_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'LS_dB_RHCP_5_zenith'}));
LS_dB_RHCP_5_zenith=instFile_struct.VALUE(indexLS_dB_RHCP_5_zenith);

indexmG_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'mG_LHCP_5'}));
mG_LHCP_5 = instFile_struct.VALUE(indexmG_LHCP_5);
indexmG_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'mG_RHCP_5'}));
mG_RHCP_5 = instFile_struct.VALUE(indexmG_RHCP_5);
indexmG_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'mG_RHCP_5_zenith'}));
mG_RHCP_5_zenith = instFile_struct.VALUE(indexmG_RHCP_5_zenith);

indexmNF_LHCP_5 = find(strcmpi(instFile_struct.RowNames,{'mNF_LHCP_5'}));
mNF_LHCP_5=instFile_struct.VALUE(indexmNF_LHCP_5);
indexmNF_RHCP_5 = find(strcmpi(instFile_struct.RowNames,{'mNF_RHCP_5'}));
mNF_RHCP_5=instFile_struct.VALUE(indexmNF_RHCP_5);
indexmNF_RHCP_5_zenith = find(strcmpi(instFile_struct.RowNames,{'mNF_RHCP_5_zenith'}));
mNF_RHCP_5_zenith=instFile_struct.VALUE(indexmNF_RHCP_5_zenith);

instrumentConfParam.F_dB_LHCP = [F_dB_LHCP_1 F_dB_LHCP_5];
instrumentConfParam.F_dB_RHCP = [F_dB_RHCP_1 F_dB_RHCP_5];
instrumentConfParam.gG_LHCP = [gG_LHCP_1 gG_LHCP_5];
instrumentConfParam.gG_RHCP = [gG_RHCP_1 gG_RHCP_5];
instrumentConfParam.gG_RHCP_zenith = [gG_RHCP_1_zenith gG_RHCP_5_zenith];
instrumentConfParam.gNF_LHCP = [gNF_LHCP_1 gNF_LHCP_5];
instrumentConfParam.gNF_RHCP = [gNF_RHCP_1 gNF_RHCP_5];
instrumentConfParam.gNF_RHCP_zenith = [gNF_RHCP_1_zenith gNF_RHCP_5_zenith];
instrumentConfParam.Gsys_dB_LHCP = [Gsys_dB_LHCP_1 Gsys_dB_LHCP_5];
instrumentConfParam.Gsys_dB_RHCP = [Gsys_dB_RHCP_1 Gsys_dB_RHCP_5];
instrumentConfParam.Gsys_dB_RHCP_zenith = [Gsys_dB_RHCP_1_zenith Gsys_dB_RHCP_5_zenith];
instrumentConfParam.LN_dB_LHCP = [LN_dB_LHCP_1 LN_dB_LHCP_5];
instrumentConfParam.LN_dB_RHCP = [LN_dB_RHCP_1 LN_dB_RHCP_5];
instrumentConfParam.LN_dB_RHCP_zenith = [LN_dB_RHCP_1_zenith LN_dB_RHCP_5_zenith];
instrumentConfParam.LS_dB_LHCP = [LS_dB_LHCP_1 LS_dB_LHCP_5];
instrumentConfParam.LS_dB_RHCP = [LS_dB_RHCP_1 LS_dB_RHCP_5];
instrumentConfParam.LS_dB_RHCP_zenith = [LS_dB_RHCP_1_zenith LS_dB_RHCP_5_zenith];
instrumentConfParam.mG_LHCP = [mG_LHCP_1 mG_LHCP_5];
instrumentConfParam.mG_RHCP = [mG_RHCP_1 mG_RHCP_5];
instrumentConfParam.mG_RHCP_zenith = [mG_RHCP_1_zenith mG_RHCP_5_zenith];
instrumentConfParam.mNF_LHCP = [mNF_LHCP_1 mNF_LHCP_5];
instrumentConfParam.mNF_RHCP = [mNF_RHCP_1 mNF_RHCP_5];
instrumentConfParam.mNF_RHCP_zenith = [mNF_RHCP_1_zenith mNF_RHCP_5_zenith];



end