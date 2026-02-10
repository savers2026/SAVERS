%% ************************************************************************
% MODULE NAME:      readDEM_5km.m
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
% The module opens a 5km DEM in geotiff format or .mat format
% and reads only a specific area of the image around the SP and defined
% by the area size given as input. The DEM must be referred to WGS84.
% The module then returns the elevation at the SP coordinates (SP was found
% on the ellipsoid).
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function surfaceElevationAtSPellipsoid = readDEM_5km(DEMinfo,coord_SP_lla_onEllipsoid)


if contains(DEMinfo.filenameDEM, 'tif')

    %read the surface elevation data
    [surfaceElevation, spatialRefInfo]=readgeoraster(DEMinfo.filenameDEM5km);
    %latitude delta of DEM
    deltaLatDEM = spatialRefInfo.CellExtentInLatitude;
    %longitude delta of DEM
    deltaLonDEM = spatialRefInfo.CellExtentInLongitude;

    %coordinates of the upper left corner of the upper left pixel
    maxLatOrigDEM=spatialRefInfo.LatitudeLimits(2);
    minLonOrigDEM=spatialRefInfo.LongitudeLimits(1);
    %coordinates of the bottom right corner of the bottom right pixel
    minLatOrigDEM=spatialRefInfo.LatitudeLimits(1);
    maxLonOrigDEM=spatialRefInfo.LongitudeLimits(2);

    %array of lat and lon of original DEM
    %these are the coordinates of the center of the pixel
    latDEM=(minLatOrigDEM+deltaLatDEM/2):deltaLatDEM:(maxLatOrigDEM-deltaLatDEM/2);
    lonDEM=(minLonOrigDEM+deltaLonDEM/2):deltaLonDEM:(maxLonOrigDEM-deltaLonDEM/2);
    latDEM=flipud(latDEM');

    %compute matrix of LAT
    DEM_lla(:,:,1)=repmat(latDEM,1,length(lonDEM));
    %compute matrix of LON
    DEM_lla(:,:,2)=repmat(lonDEM,length(latDEM),1);
    DEM_lla(:,:,3) = surfaceElevation;
    clear surfaceElevation

else
    load(DEMinfo.filenameDEM)
    deltaLatDEM=abs(DEM_lla(1,1,1)-DEM_lla(2,1,1));
    deltaLonDEM=abs(DEM_lla(1,1,2)-DEM_lla(1,2,2));
end

%look for the DEM elevation corresponding at the SP coordinates found on the ellipsoid
indexSPlatInDEM = find((abs(DEM_lla(:,1,1) - coord_SP_lla_onEllipsoid(1))) < deltaLatDEM/2);
indexSPlonInDEM = find((abs(DEM_lla(1,:,2) - coord_SP_lla_onEllipsoid(2))) < deltaLonDEM/2);
surfaceElevationAtSPellipsoid = double(DEM_lla(indexSPlatInDEM,indexSPlonInDEM,3));

return
