%% ************************************************************************
% MODULE NAME:      biogeoparamAtSP.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      September 2022 by Laura Dente
% Main updates:          
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% This module returns in a structure the bio geo parameters used for
% for the simulation at the nominal specular point: cover class (after
% interpolation and rotation over the simulation grid), soil status
% (normal/flooded/frozen), soil moisture and soil roughness, LAi and
% biomass.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function bioGeoParametersAtSP=biogeoparamAtSP(infoAtSPforRefFile,coverMap_IDs,terrainParameters,vegetationParameters)

% %look for the cover class where the SP falls
% deltaLatMap=abs(landCoverMap(1,1,1)-landCoverMap(2,1,1));
% deltaLonMap=abs(landCoverMap(1,1,2)-landCoverMap(1,2,2));
numSP=size(infoAtSPforRefFile.coverClassAtSP,2);
for i=1:numSP
%     indexSPlatInMap = find((abs(landCoverMap(:,1,1)-SPlatlonAlt(i,1)))<deltaLatMap/2);
%     indexSPlonInMap = find((abs(landCoverMap(1,:,2)-SPlatlonAlt(i,2)))<deltaLonMap/2);
    classOfSPinMap = infoAtSPforRefFile.coverClassAtSP(i);
    %look for the bio/geo parameters given as input for the SP class
    indexClass=find(coverMap_IDs==classOfSPinMap);
    %Amir direct signal - should be not commented

    soilStatusAtSP(i)=terrainParameters(indexClass,1);
    soilRoughnessAtSP(i)=terrainParameters(indexClass,2);
    soilMoistureAtSP(i)=terrainParameters(indexClass,5);
    LAIatSP(i)=vegetationParameters(indexClass,3);
    BioAtSP(i)=vegetationParameters(indexClass,7);

%     soilStatusAtSP(i)=0;
%     soilRoughnessAtSP(i)=0;
%     soilMoistureAtSP(i)=0;
%     LAIatSP(i)=0;
%     BioAtSP(i)=0;

end

bioGeoParametersAtSP.classOfSPinCoverMap=infoAtSPforRefFile.coverClassAtSP;
bioGeoParametersAtSP.soilStatusAtSP=soilStatusAtSP;
bioGeoParametersAtSP.soilRoughnessAtSP=soilRoughnessAtSP;
bioGeoParametersAtSP.soilMoistureAtSP=soilMoistureAtSP;
bioGeoParametersAtSP.LAIatSP=LAIatSP;
bioGeoParametersAtSP.BioAtSP=BioAtSP;

end
