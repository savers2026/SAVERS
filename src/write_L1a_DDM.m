%% ************************************************************************
% MODULE NAME:      write_L1a_DDM.m
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
% The module writes the DDMs.nc L1a product of the Hydrognss E2ES.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function write_L1a_DDM(TrackID,geoSYSp,IntegrationMidPointTime,DDM_LR1_linear,DDM_LR5_linear,DDM_RR1_linear,DDM_RR5_linear,...
                                    outpathname,foldername, CoherentMeasurementI_LR1, CoherentMeasurementI_LR5, CoherentMeasurementQ_LR1, CoherentMeasurementQ_LR5, IntegrationMidPointTime_coe,TrackID_cohChannel)



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

Delay=(0:geoSYSp.nx_Rdelay-1)';
% Delay_linear=10.^(Delay/10);
[a,b]= size(Delay);
Doppler=(-fix(geoSYSp.ny_Dfreq/2) : round(geoSYSp.ny_Dfreq/2)-1)';

% Doppler_linear=10.^(Doppler/10);
[a1,b1]=size(Doppler);
% IntegrationMidPointTime=datenum(1980, 1, 6, 0 , 0 ,0)+GNSS_Week.*7+GNSS_SoW./86400;
[a2 b2]=size(IntegrationMidPointTime);

% DDM_linear  = outputSimulations.sigPower_Noise_LR_lin;
% DDM_linear2 = outputSimulations.sigPower_Noise_RR_lin;

maxddm_LR1  = max(max(DDM_LR1_linear));
maxddm_RR1  = max(max(DDM_RR1_linear));
maxddm_LR1  = squeeze(maxddm_LR1);
maxddm_RR1  = squeeze(maxddm_RR1);

maxddm_LR5  = max(max(DDM_LR5_linear));
maxddm_RR5  = max(max(DDM_RR5_linear));
maxddm_LR5  = squeeze(maxddm_LR5);
maxddm_RR5  = squeeze(maxddm_RR5);

%in case of simulation for L1 or E1 only, set max for L5 and E5 to 1 (to avoid problems
%in the normalization to max)
if maxddm_LR5(1) == 0
    maxddm_LR5=maxddm_LR5+1.0;
    maxddm_RR5=maxddm_LR5+1.0;
    numch=2; %no of channels
else
    numch=4; %no of channels
end

[ad bd]=size(maxddm_LR1);
for ki=1:ad
    ddm_LR1_digNum(:,:,ki)= DDM_LR1_linear(:,:,ki).*65535./(maxddm_LR1(ki));
    ddm_RR1_digNum(:,:,ki)= DDM_RR1_linear(:,:,ki).*65535./(maxddm_RR1(ki));
    ddm_LR5_digNum(:,:,ki)= DDM_LR5_linear(:,:,ki).*65535./(maxddm_LR5(ki));
    ddm_RR5_digNum(:,:,ki)= DDM_RR5_linear(:,:,ki).*65535./(maxddm_RR5(ki));
end
%   max_new=max(max(ddm)); %to verify the max value of the linear ddm

ddm_LR1_digNum = uint16(ddm_LR1_digNum);
ddm_RR1_digNum = uint16(ddm_RR1_digNum);
ddm_LR5_digNum = uint16(ddm_LR5_digNum);
ddm_RR5_digNum = uint16(ddm_RR5_digNum);
% ddm_t1=zeros(size(ddm1));
% ddm_t2=zeros(size(ddm2));
% ddm_t=[];
%for t=1:a2
ddm_LR1_digNum_transp = pagetranspose(ddm_LR1_digNum);
ddm_RR1_digNum_transp = pagetranspose(ddm_RR1_digNum);
ddm_LR5_digNum_transp = pagetranspose(ddm_LR5_digNum);
ddm_RR5_digNum_transp = pagetranspose(ddm_RR5_digNum);
%end

InstantaneousPower_LR1   = (CoherentMeasurementI_LR1.^2 + CoherentMeasurementQ_LR1.^2);
InstantaneousPower_LR5   = (CoherentMeasurementI_LR5.^2 + CoherentMeasurementQ_LR5.^2);

%creatingtheddmfile
% ncid = netcdf.create('DDMs.nc','NETCDF4');

ncid = netcdf.create(strjoin([outpathname '\' 'DDMs.nc'],''),'netcdf4');
netcdf.defDim (ncid, 'Time' , netcdf.getConstant ( 'NC_UNLIMITED' ));
date_created=datetime('now','Format','yyyy-MM-dd''T''hh:mm:ss');
date_created=string(date_created);
L1a_ProcessorVersion=single(1.0);
%global attributes
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'Name','HydroGNSS L1a DDM File');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'date_created',date_created);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'License','License: HydroGNSS GNSS-R Dataset');
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'FileIDCode',foldername);
netcdf.putAtt(ncid, netcdf.getConstant("NC_GLOBAL"),'L1a_ProcessorVersion',L1a_ProcessorVersion); % changed by ansha

for n= 1:nooftracks
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
    for m= 1:numch
        if m==1
            ddm=ddm_LR1_digNum_transp;
            CoherentMeasurementI = CoherentMeasurementI_LR1;
            CoherentMeasurementQ = CoherentMeasurementQ_LR1;
            InstantaneousPower = InstantaneousPower_LR1;
        elseif m==2
            ddm=ddm_RR1_digNum_transp;
            CoherentMeasurementI = zeros(length(IntegrationMidPointTime_coe),1);
            CoherentMeasurementQ = zeros(length(IntegrationMidPointTime_coe),1);
            InstantaneousPower = zeros(length(IntegrationMidPointTime_coe),1);
        elseif m==3
            ddm=ddm_LR5_digNum_transp;
            CoherentMeasurementI = CoherentMeasurementI_LR5;
            CoherentMeasurementQ = CoherentMeasurementQ_LR5;
            InstantaneousPower = InstantaneousPower_LR5;
        else
            ddm=ddm_RR5_digNum_transp;
            CoherentMeasurementI = zeros(length(IntegrationMidPointTime_coe),1);
            CoherentMeasurementQ = zeros(length(IntegrationMidPointTime_coe),1);
            InstantaneousPower = zeros(length(IntegrationMidPointTime_coe),1);
        end
        trackid_ch(m) = netcdf.defGrp(trackid(n),['Channel' num2str(m)]);
        netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'ChannelID',uint16(m));
        netcdf.putAtt(trackid_ch(m), netcdf.getConstant("NC_GLOBAL"),'TrackID',uint32(n-1));
        trackid_n1(m) = netcdf.defGrp(trackid_ch(m),'Incoherent');
        dimid = netcdf.defDim(trackid_n1(m),'Delay',a);
        var1=netcdf.defVar(trackid_n1(m),'Delay','NC_DOUBLE',dimid);
        netcdf.putVar(trackid_n1(m),var1,Delay);
        dimid1 = netcdf.defDim(trackid_n1(m),'Doppler',a1);
        var2=netcdf.defVar(trackid_n1(m),'Doppler','NC_DOUBLE',dimid1);
        netcdf.putVar(trackid_n1(m),var2,Doppler);
        index = find(newtrackid(n)==TrackID);
        [l1 l2]=size(index);
        dimid2=netcdf.defDim(trackid_n1(m),'DateTimeTrack',l1);
        var3=netcdf.defVar(trackid_n1(m),'IntegrationMidPointTime','NC_DOUBLE',dimid2);
        netcdf.putVar(trackid_n1(m),var3,IntegrationMidPointTime(index));
        dimid4=([dimid dimid1 dimid2]);
        var4=netcdf.defVar(trackid_n1(m),'DDM','NC_USHORT',dimid4);
        netcdf.putVar(trackid_n1(m),var4,ddm(:,:,index));

        trackid_n2(m) = netcdf.defGrp(trackid_ch(m),'Coherent');
        indexCohChannel = find(newtrackid(n)==TrackID_cohChannel);
        l1coe = length(indexCohChannel);

        fact_coe = 1;

        dimid3 = netcdf.defDim(trackid_n2(m),'CoherentDateTimeTrack',l1coe);
        var5   = netcdf.defVar(trackid_n2(m),'CoherentIntegrationMidPointTime','NC_DOUBLE',dimid3);
        netcdf.putVar(trackid_n2(m), var5, IntegrationMidPointTime_coe(indexCohChannel));
        var6   = netcdf.defVar(trackid_n2(m),'CoherentMeasurementI','NC_INT',dimid3);
        netcdf.putVar(trackid_n2(m),var6,fact_coe.*CoherentMeasurementI(indexCohChannel));
        var7   = netcdf.defVar(trackid_n2(m),'CoherentMeasurementQ','NC_INT',dimid3);
        netcdf.putVar(trackid_n2(m),var7,fact_coe.*CoherentMeasurementQ(indexCohChannel));
        var8   = netcdf.defVar(trackid_n2(m),'InstantaneousPower','NC_UINT',dimid3);
        netcdf.putVar(trackid_n2(m),var8,fact_coe.*InstantaneousPower(indexCohChannel));

    end
end
netcdf.close(ncid);
% DDM_HydroGNSS=ncinfo('DDMs.nc');
end

