%% ************************************************************************
% MODULE NAME:      WriteSixHours.m
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
% Segmentation of the simulator output in six-hours segments/directories.
% This module calls all the routines to write the L1a products.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [UTCDateMidPoint] = WriteSixHours(satellite,TRAp,Gsys_dB,GNSS_Week,GNSS_SoW,Rx,VRx,Tx,VTx,SPlatlonAlt,PRN,TrackID,TrackID_cohChannel,SVN,outputSimulations,name_data_tag,geoSYSp,outfile_directory)
 
  IntegrationMidPointTime = datenum(1980, 1, 6, 0 , 0 ,0)+(GNSS_Week).*7+GNSS_SoW./86400 ;
  % % Identify number of specular point in the whole set of tracks 
  [ai,bi]=size(IntegrationMidPointTime);

if ai==1
    IntegrationMidPointTime=IntegrationMidPointTime';
    PRN=PRN';
    TrackID=uint32(TrackID)';
    SVN=SVN';
    GNSS_SoW=GNSS_SoW';
    GNSS_Week=GNSS_Week';
end
% IntegrationMidPointTime=IntegrationMidPointTime(1:21000) ; 
 
% [SizeL1A b]=size(IntegrationMidPointTime) ; 

% YearMidPoint=[] ; 
% MonthMidPoint=[] ;
% DayMidPoint=[] ; 
% HourMidPoint=[]  ; 
% Here I am assuming you have the data in a vector with the same size of 
% IntegrationMidPoint time. I am filling it with dummy values =0 
% metadata=zeros(SizeL1A, 1) ; 
% Identify year, month, day and six hour block to be separated in different
% directories
UTCDateMidPoint=datetime(IntegrationMidPointTime, 'ConvertFrom', 'datenum') ; 
[B,I]=sort(UTCDateMidPoint);
UTCDateMidPoint=UTCDateMidPoint(I) ; 
GNSS_SoW=GNSS_SoW(I);
GNSS_Week=GNSS_Week(I);
Rx=Rx(I,:);
VRx=VRx(I,:);
Tx=Tx(I,:);
VTx=VTx(I,:);

PRN=PRN(I);
TrackID=TrackID(I);
SVN=SVN(I);
%  TrackIDOrbit=TrackIDOrbit(I) ; 



IntegrationMidPointTime = IntegrationMidPointTime(I);

% DDMOutputNumericalScaling_noNoise_LR=outputSimulations.max_sigPower_noNoise_LR(I);
% DDMOutputNumericalScaling_noNoise_RR=outputSimulations.max_sigPower_noNoise_RR(I);
% DDMOutputNumericalScaling_Noise_LR=outputSimulations.max_sigPower_Noise_LR(I);
% DDMOutputNumericalScaling_Noise_RR=outputSimulations.max_sigPower_Noise_RR(I);
DDMOutputNumericalScaling_NLblackbody_DOWNmax_LR_1 = outputSimulations.NLblackbody_DOWNmax_LR_1(I);
DDMOutputNumericalScaling_NLblackbody_DOWNmax_LR_5 = outputSimulations.NLblackbody_DOWNmax_LR_5(I);
DDMOutputNumericalScaling_NLblackbody_DOWNmax_RR_1 = outputSimulations.NLblackbody_DOWNmax_RR_1(I);
DDMOutputNumericalScaling_NLblackbody_DOWNmax_RR_5 = outputSimulations.NLblackbody_DOWNmax_RR_5(I);
DDMOutputNumericalScaling_NLblackbody_UPmax_1      = outputSimulations.NLblackbody_UPmax_1(I);
% DDMOutputNumericalScaling_NLblackbody_UPmax_5      = outputSimulations.NLblackbody_UPmax_5(I);

SignalPower_ds_1   = outputSimulations.RR1_UP_Power_Fluc(I);
%SignalPower_ds_5  = outputSimulations.RR5_UP_Power_Fluc(I);
NoisePower_ds_1    = outputSimulations.PN_star_1(I);
NoisePower_ds_5    = outputSimulations.PN_star_5(I);
BlackBodyPower_ds  = zeros(size(SignalPower_ds_1)); 
DDM_NLBB_down_LR_1 = outputSimulations.NLblackbody_DOWN_LR_1(:,:,I);
DDM_NLBB_down_LR_5 = outputSimulations.NLblackbody_DOWN_LR_5(:,:,I);
DDM_NLBB_down_RR_1 = outputSimulations.NLblackbody_DOWN_RR_1(:,:,I);
DDM_NLBB_down_RR_5 = outputSimulations.NLblackbody_DOWN_RR_5(:,:,I);
DDM_NLBB_up_1      = outputSimulations.NLblackbody_UP_1(:,:,I);
DDM_NLBB_up_5      = outputSimulations.NLblackbody_UP_5(:,:,I);

DDMPixelValueNoiseWholeDDM_down_LR_1 = outputSimulations.NLblackbody_DOWNmean_LR_1(I);
DDMPixelValueNoiseWholeDDM_down_LR_5 = outputSimulations.NLblackbody_DOWNmean_LR_5(I);
DDMPixelValueNoiseWholeDDM_down_RR_1 = outputSimulations.NLblackbody_DOWNmean_RR_1(I);
DDMPixelValueNoiseWholeDDM_down_RR_5 = outputSimulations.NLblackbody_DOWNmean_RR_5(I);
DDMPixelValueNoiseWholeDDM_up_1      = outputSimulations.NLblackbody_UPmean_1(I);
DDMPixelValueNoiseWholeDDM_up_5      = outputSimulations.NLblackbody_UPmean_5(I);

max_sigPower_Noise_LR1 = outputSimulations.max_sigPower_Noise_LR1(I);
max_sigPower_Noise_LR5 = outputSimulations.max_sigPower_Noise_LR5(I);
max_sigPower_Noise_RR1 = outputSimulations.max_sigPower_Noise_RR1(I);
max_sigPower_Noise_RR5 = outputSimulations.max_sigPower_Noise_RR5(I);
%%variables assigning from structures
%modified by Laura
% ReflectionHeight=geoSYSp.SPAlt_series(I);
% SatSpecularPointLat=geoSYSp.SPlat_series(I);
% SatSpecularPointLon=geoSYSp.SPlon_series(I);
ReflectionHeight    = SPlatlonAlt(I,3);
SatSpecularPointLat = SPlatlonAlt(I,1);
SatSpecularPointLon = SPlatlonAlt(I,2);

DDM_LR1_linear = outputSimulations.sigPower_Noise_LR1_lin(:,:,I);
DDM_LR5_linear = outputSimulations.sigPower_Noise_LR5_lin(:,:,I);
DDM_RR1_linear = outputSimulations.sigPower_Noise_RR1_lin(:,:,I);
DDM_RR5_linear = outputSimulations.sigPower_Noise_RR5_lin(:,:,I);

OnBoardDelayReferencedToDirect_1 = outputSimulations.delayAtSPReferencedToDirect_1_microsec(I)*10^3; %nanosec
OnBoardDelayReferencedToDirect_5 = outputSimulations.delayAtSPReferencedToDirect_5_microsec(I)*10^3; %nanosec


[SizeL1A,r]   = size(UTCDateMidPoint);
YearMidPoint  = year(UTCDateMidPoint) ; 
MonthMidPoint = month(UTCDateMidPoint) ;
DayMidPoint   = day(UTCDateMidPoint) ; 
HourMidPoint  = hour(UTCDateMidPoint)  ; 
% Create the Hour name as 00, 06, 12, or 18
HourMidPoint  = 6*fix(HourMidPoint/6) ; 
% Create name of directory for each specular point 
Path_L1Adata  = [] ; 
formatSpec = '%02u' ; 

for ii=1:SizeL1A
    read=[num2str(YearMidPoint(ii)), '-', num2str(MonthMidPoint(ii),formatSpec),...
        '/', num2str(DayMidPoint(ii), formatSpec), '/H',num2str(HourMidPoint(ii), formatSpec) ] ;
    Path_L1Adata=char(Path_L1Adata, read) ;
end

Path_L1Adata = string(Path_L1Adata) ;
Path_L1Adata = Path_L1Adata(2:end) ; 
% Identify blocks of specular point that should be filled in a unique 
% direcotry ad they pertain to the same six hour block
[Directories,a,b] = unique(Path_L1Adata);
[NumbSixHours,c]  = size(Directories) ;
% Segment the whole file in input into six hours block by identifying the
% index (variable IndexSixHours_ii) of point pertaining to each six hour
% block (block ii of NumSixHours). Then write the subset file into the 
% proper directory
% TrackIDorbit=[] ; % Added by Mauro to create Look Up Table with TrackID in metadta and from Orbit file

%%
 CoherentMeasurementI_1 = real(outputSimulations.Atotal_1);
 CoherentMeasurementI_5 = real(outputSimulations.Atotal_5);
 CoherentMeasurementQ_1 = imag(outputSimulations.Atotal_1);
 CoherentMeasurementQ_5 = imag(outputSimulations.Atotal_5);
 if CoherentMeasurementI_1 == 0
     IntegrationMidPointTime_coe = zeros(1,length(CoherentMeasurementI_1));
 else
    M                           = geoSYSp.incoIntTime./geoSYSp.cohIntTime;
    IntegrationMidPointTime_coe = linspace(IntegrationMidPointTime(1) - geoSYSp.incoIntTime/(24*60*60)/1000/2, IntegrationMidPointTime(end) + geoSYSp.incoIntTime/(24*60*60)/1000/2, M*length(IntegrationMidPointTime));  
 end
%%

for ii=1:NumbSixHours 
   disp(['Six hour block  ', num2str(ii), ' of ', num2str(NumbSixHours)]) ;
    IndexSixHours_ii= find(b == ii) ; 

    foldername  = append(Directories(ii)+'/');
    outpathname = ([outfile_directory '\' name_data_tag '\' satellite '\' 'DataRelease\L1A_L1B\' convertStringsToChars(Directories(ii))]);
    outpathname = replace(outpathname,'/','\');
    mkdir(strjoin(outpathname,''));
%     selectAndReadGNSSRdata_metadata(delayOffset,GNSS_SoW(IndexSixHours_ii),GNSS_Week(IndexSixHours_ii),...
%                           Rx(IndexSixHours_ii,:),VRx(IndexSixHours_ii,:),Tx(IndexSixHours_ii,:),VTx(IndexSixHours_ii,:),...
%                           PRN(IndexSixHours_ii),TrackID(IndexSixHours_ii),TrackIDOrbit(IndexSixHours_ii),SVN(IndexSixHours_ii),...
%                           name_data_tag,geoSYSp,IntegrationMidPointTime(IndexSixHours_ii),ReflectionHeight(IndexSixHours_ii),...
%                           SatSpecularPointLat(IndexSixHours_ii),SatSpecularPointLon(IndexSixHours_ii),...
%                           max_sigPower_Noise_LR(IndexSixHours_ii),max_sigPower_Noise_RR(IndexSixHours_ii),outpathname,foldername);
   write_L1a_metadata(satellite,TRAp,GNSS_SoW(IndexSixHours_ii),GNSS_Week(IndexSixHours_ii),...
                            Rx(IndexSixHours_ii,:),VRx(IndexSixHours_ii,:),Tx(IndexSixHours_ii,:),VTx(IndexSixHours_ii,:),...
                            PRN(IndexSixHours_ii),TrackID(IndexSixHours_ii),SVN(IndexSixHours_ii),name_data_tag,geoSYSp,IntegrationMidPointTime(IndexSixHours_ii),...
                            ReflectionHeight(IndexSixHours_ii),SatSpecularPointLat(IndexSixHours_ii),SatSpecularPointLon(IndexSixHours_ii),...
                            OnBoardDelayReferencedToDirect_1(IndexSixHours_ii),OnBoardDelayReferencedToDirect_5(IndexSixHours_ii),...
                            max_sigPower_Noise_LR1(IndexSixHours_ii),max_sigPower_Noise_LR5(IndexSixHours_ii),max_sigPower_Noise_RR1(IndexSixHours_ii),max_sigPower_Noise_RR5(IndexSixHours_ii),...
                            outpathname,foldername, IntegrationMidPointTime_coe);
   write_L1a_DDM(TrackID(IndexSixHours_ii),geoSYSp,IntegrationMidPointTime(IndexSixHours_ii),...
                            DDM_LR1_linear(:,:,IndexSixHours_ii),DDM_LR5_linear(:,:,IndexSixHours_ii),DDM_RR1_linear(:,:,IndexSixHours_ii),DDM_RR5_linear(:,:,IndexSixHours_ii),...
                            outpathname,foldername, CoherentMeasurementI_1,CoherentMeasurementI_5, CoherentMeasurementQ_1,CoherentMeasurementQ_5,...
                            IntegrationMidPointTime_coe,TrackID_cohChannel);
   write_L1a_bbZenith(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_UPmax_1(IndexSixHours_ii), DDMPixelValueNoiseWholeDDM_up_1(IndexSixHours_ii),...
                            DDM_NLBB_up_1(:,:,IndexSixHours_ii), PRN(IndexSixHours_ii), outpathname, foldername);
   
   if geoSYSp.indexTxSignal == 0 || geoSYSp.indexTxSignal == 2
         write_L1a_bbNadirLR1(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_LR_1(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_LR_1(IndexSixHours_ii),...
                            DDM_NLBB_down_LR_1(:,:,IndexSixHours_ii), PRN(IndexSixHours_ii), outpathname, foldername);
         write_L1a_bbNadirRR1(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_RR_1(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_RR_1(IndexSixHours_ii),...
                            DDM_NLBB_down_RR_1(:,:,IndexSixHours_ii),PRN(IndexSixHours_ii), outpathname, foldername);
   else
         write_L1a_bbNadirLR1(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_LR_1(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_LR_1(IndexSixHours_ii),...
                            DDM_NLBB_down_LR_1(:,:,IndexSixHours_ii), PRN(IndexSixHours_ii), outpathname, foldername);
         write_L1a_bbNadirRR1(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_RR_1(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_RR_1(IndexSixHours_ii),...
                            DDM_NLBB_down_RR_1(:,:,IndexSixHours_ii),PRN(IndexSixHours_ii), outpathname, foldername);
         write_L1a_bbNadirLR5(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_LR_5(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_LR_5(IndexSixHours_ii),...
                            DDM_NLBB_down_LR_5(:,:,IndexSixHours_ii), PRN(IndexSixHours_ii), outpathname, foldername);
        write_L1a_bbNadirRR5(IntegrationMidPointTime(IndexSixHours_ii), geoSYSp, Rx(IndexSixHours_ii,:),...
                            DDMOutputNumericalScaling_NLblackbody_DOWNmax_RR_5(IndexSixHours_ii),...
                            DDMPixelValueNoiseWholeDDM_down_RR_5(IndexSixHours_ii),...
                            DDM_NLBB_down_RR_5(:,:,IndexSixHours_ii), PRN(IndexSixHours_ii), outpathname, foldername);
  end
 



%  write_L1a_bbNadir(IntegrationMidPointTime(IndexSixHours_ii),geoSYSp,Rx(IndexSixHours_ii,:),...
%                             DDMOutputNumericalScaling_NLblackbody_DOWNmax_1(IndexSixHours_ii),DDMOutputNumericalScaling_NLblackbody_DOWNmax_5(IndexSixHours_ii), ...
%                             DDMPixelValueNoiseWholeDDM_down_1(IndexSixHours_ii),DDMPixelValueNoiseWholeDDM_down_5(IndexSixHours_ii),...
%                             DDM_NLBB_down_1(:,:,IndexSixHours_ii),DDM_NLBB_down_5(:,:,IndexSixHours_ii),PRN(IndexSixHours_ii),outpathname,foldername);
   write_L1a_DSpow(TRAp,Gsys_dB,Rx(IndexSixHours_ii,:),VRx(IndexSixHours_ii,:),Tx(IndexSixHours_ii,:),VTx(IndexSixHours_ii,:),...
                            SignalPower_ds_1(IndexSixHours_ii),NoisePower_ds_1(IndexSixHours_ii),BlackBodyPower_ds(IndexSixHours_ii),...
                            PRN(IndexSixHours_ii),SVN(IndexSixHours_ii),IntegrationMidPointTime(IndexSixHours_ii),geoSYSp,outpathname)
end

end