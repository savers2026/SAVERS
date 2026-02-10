%% ************************************************************************
% MODULE NAME:      siteDependentInputs.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome
% Original version:      Sep 2019 by Laura Dente
% Main updates:          
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module returns the DEm filenames and the area limits according to the
% scenario and site selcted by the user.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [filenameDEM30forPlot,filenameDEM,filenameDEM_flooded,SPcoodinatesLimits,scenario]=siteDependentInputs(soilmoisture,forest,freeze_thaw,wetlands,simResolInX,simResolInY,configParam)

%associate coordinates limits and DEM corresponding to the site
switch soilmoisture
    %Tibesti mountains
    case 'Tibesti Mountains'
        filenameDEM3arcsec=configParam.Chad.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Chad.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Chad.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Chad.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Chad.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Chad.DEM30arcsec;% 900 m resolution
        %areaDEM=[16.0 25.0 14.0 22.0];
        SPcoodinatesLimits=configParam.Chad.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=1;        

    %San Luis valley
    case 'San Luis Valley'
        filenameDEM3arcsec=configParam.SLV.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.SLV.DEM9arcsec;
        filenameDEM30arcsec=configParam.SLV.DEM30arcsec;
        filenameDEM3arcsec_flooded=configParam.SLV.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.SLV.DEM9arcsec;
        filenameDEM30arcsec_flooded=configParam.SLV.DEM30arcsec;        
        %areaDEM=[33.5 41.5 -110.0 -102.0];
        SPcoodinatesLimits=configParam.SLV.SPcoordLimits; %minLat maxLat minLon maxLo
        scenario=1;
    case 'Site 3'
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Input error: This site selection is not active yet.';...
                'The simulation will be stopped.';...
                'Please run again HSAVERS selecting a different site'},'INPUT ERROR MESSAGE',CreateStruct))
        close all
        return        
end

switch forest
    %Gabon
    case 'Gabon'
        filenameDEM3arcsec=configParam.Gabon.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Gabon.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Gabon.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Gabon.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Gabon.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Gabon.DEM30arcsec;% 900 m resolution
        %areaDEM=[-6.0 2.0 9.0 17.0];
        SPcoodinatesLimits=configParam.Gabon.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=2;

    % Sodankyla and Saariselka
    case 'Finland'
        filenameDEM3arcsec=configParam.Finland.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Finland.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Finland.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Finland.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Finland.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Finland.DEM30arcsec;% 900 m resolution
        %areaDEM=[63.0 72.0 22.0 32.0];
        SPcoodinatesLimits=configParam.Finland.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=2;

    case 'Site 3'
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Input error: This site selection is not active yet.';...
                'The simulation will be stopped.';...
                'Please run again HSAVERS selecting a different site'},'INPUT ERROR MESSAGE',CreateStruct))
        close all
        return        

end

switch freeze_thaw

    % Sodankyla and Saariselka
    case 'Finland'
        filenameDEM3arcsec=configParam.Finland.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Finland.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Finland.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Finland.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Finland.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Finland.DEM30arcsec;% 900 m resolution
        %areaDEM=[63.0 72.0 22.0 32.0];
        SPcoodinatesLimits=configParam.Finland.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=3;

    % Chersky, Siberia
    case 'Chersky'
        filenameDEM3arcsec=configParam.Chersky.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Chersky.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Chersky.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Chersky.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Chersky.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Chersky.DEM30arcsec;% 900 m resolution
        %areaDEM=[64.0 72.0 158.0 166.0];
        SPcoodinatesLimits=configParam.Chersky.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=3;

    case 'Site 3'
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Input error: This site selection is not active yet.';...
                'The simulation will be stopped.';...
                'Please run again HSAVERS selecting a different site'},'INPUT ERROR MESSAGE',CreateStruct))
        close all
        return        

end

switch wetlands
    %Yucatan lake
    case 'Yucatan Lake'
        filenameDEM3arcsec=configParam.Yucatan.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Yucatan.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Yucatan.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Yucatan.DEM3arcsec_flooded;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Yucatan.DEM9arcsec_flooded;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Yucatan.DEM30arcsec_flooded;% 900 m resolution
        %areaDEM=[28.0 36.0 -95.0 -87.0];
        SPcoodinatesLimits=configParam.Yucatan.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=4;        

    case 'Everglades'
        filenameDEM3arcsec=configParam.Everglades.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec=configParam.Everglades.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec=configParam.Everglades.DEM30arcsec;% 900 m resolution
        filenameDEM3arcsec_flooded=configParam.Everglades.DEM3arcsec;% 90 m resolution
        filenameDEM9arcsec_flooded=configParam.Everglades.DEM9arcsec;% 290 m resolution
        filenameDEM30arcsec_flooded=configParam.Everglades.DEM30arcsec;% 900 m resolution
        %areaDEM=[22.0 30.0 -85.0 -77.0];
        SPcoodinatesLimits=configParam.Everglades.SPcoordLimits; %minLat maxLat minLon maxLon
        scenario=4;        
        
    case 'Site 3'
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        uiwait(msgbox({'Input error: This site selection is not active yet.';...
                'The simulation will be stopped.';...
                'Please run again HSAVERS selecting a different site'},'INPUT ERROR MESSAGE',CreateStruct))
        close all
        return

end

if simResolInX<290 && simResolInY<290
    filenameDEM=filenameDEM3arcsec;
    filenameDEM_flooded=filenameDEM3arcsec_flooded;
elseif simResolInX<900 && simResolInX>=290 && simResolInY<900 && simResolInY>=290
    filenameDEM=filenameDEM9arcsec;
    filenameDEM_flooded=filenameDEM9arcsec_flooded;
elseif simResolInX>=900 && simResolInY>=900
    filenameDEM=filenameDEM30arcsec;
    filenameDEM_flooded=filenameDEM30arcsec_flooded;
else
    errordlg('Error in the input simulation resolution','Input error');
    return
end

filenameDEM30forPlot=filenameDEM30arcsec;

end