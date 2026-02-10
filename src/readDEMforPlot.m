%% ************************************************************************
% MODULE NAME:      readDEMforPlot.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome
% Original version:      2019 by Laura Dente
% Main updates:          2022 by L. Dente
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module opens a DEM in geotiff format or .mat format (obtained by a mosaic
% of SRTM DEM images) and reads only a specific area defined
% by the area size given as input. The DEM must be referred to WGS84.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [DEM_lla_resized,deltaLatDEM,deltaLonDEM] = readDEMforPlot(filenameDEM,areaSize)

if contains(filenameDEM, 'tif')

    %read the surface elevation data
    [surfaceElevation_zero, spatialRefInfo] = readgeoraster(filenameDEM);

    %latitude delta of DEM
    deltaLatDEM  =spatialRefInfo.CellExtentInLatitude;
    %longitude delta of DEM
    deltaLonDEM = spatialRefInfo.CellExtentInLongitude;

    %coordinates of the upper left corner of the upper left pixel
    maxLat=spatialRefInfo.LatitudeLimits(2);
    minLon=spatialRefInfo.LongitudeLimits(1);
    %coordinates of the bottom right corner of the bottom right pixel
    minLat=spatialRefInfo.LatitudeLimits(1);
    maxLon=spatialRefInfo.LongitudeLimits(2);

    %DEM size
    rasterWidth=spatialRefInfo.RasterSize(2);
    rasterHeight=spatialRefInfo.RasterSize(1);
    %array of lat and lon of original DEM
    %these are the coordinates of the center of the pixel
    latDEM_zero=(minLat+deltaLatDEM/2):deltaLatDEM:(maxLat-deltaLatDEM/2);
    lonDEM_zero=(minLon+deltaLonDEM/2):deltaLonDEM:(maxLon-deltaLonDEM/2);
    %when ENVI makes the mosaic, it introduces a column of zeros on the right of
    %the image and a line of zeros at the bottom of the image.
    %These column and line are here removed from the lat and lon array
    latDEM=latDEM_zero(2:rasterHeight);
    lonDEM=lonDEM_zero(1:rasterWidth-1);
    latDEM=flipud(latDEM');

    clear latDEM_zero
    clear lonDEM_zero

    %the column and line of zeros introduced by ENVI are removed
    surfaceElevation=surfaceElevation_zero(1:rasterHeight-1,1:rasterWidth-1);
    clear surfaceElevation_zero

    %compute matrix of LAT
    DEM_lla(:,:,1)=repmat(latDEM,1,length(lonDEM));
    %compute matrix of LON
    DEM_lla(:,:,2)=repmat(lonDEM,length(latDEM),1);
    DEM_lla(:,:,3) = surfaceElevation;
    clear surfaceElevation

else
    load(filenameDEM)
    deltaLatDEM=abs(DEM_lla(1,1,1)-DEM_lla(2,1,1));
    deltaLonDEM=abs(DEM_lla(1,1,2)-DEM_lla(1,2,2));
end
    
%resize the lan and lon array according the size of the selected area
indexLatDEMresized= find(DEM_lla(:,1,1)>=areaSize(1) & DEM_lla(:,1,1)<=areaSize(2));
indexLonDEMresized=find(DEM_lla(1,:,2)>=areaSize(3) & DEM_lla(1,:,2)<=areaSize(4));
DEM_lla_resized(:,:,1)=DEM_lla(indexLatDEMresized,indexLonDEMresized,1);
DEM_lla_resized(:,:,2)=DEM_lla(indexLatDEMresized,indexLonDEMresized,2);
DEM_lla_resized(:,:,3)=DEM_lla(indexLatDEMresized,indexLonDEMresized,3);
clear DEM_lla
clear indexLatDEMresized
clear indexLonDEMresized

return
