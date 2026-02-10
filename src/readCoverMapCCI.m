%% ************************************************************************
% MODULE NAME:      readCoverMapCCI.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome
% Original version:      Oct 2021 by Laura Dente
% Main updates:          2022 by L. Dente
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module opens the ESA CCI land cover map, in netcdf format,
% reads only a specific area of the image defined
% by the input SP coordinates.
% The area should cover the integration area of all the SPs.
% The land cover is simplified to include only a limited number of classes.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [landCoverMap,dominantClass_ID,change_dominant_to_land_flag]=readCoverMapCCI(coverMapFilename, SP_lla_onEllipsoid,sizeInputMapWindow)

%define the area to be read in the file
maxLatArea=max(SP_lla_onEllipsoid(:,1))+sizeInputMapWindow;
minLatArea=min(SP_lla_onEllipsoid(:,1))-sizeInputMapWindow;
maxLonArea=max(SP_lla_onEllipsoid(:,2))+sizeInputMapWindow;
minLonArea=min(SP_lla_onEllipsoid(:,2))-sizeInputMapWindow;
%look for the indexes in the file,, corresponding to the area to be read
latLandCoverMap=ncread(coverMapFilename, 'lat');
lonLandCoverMap=ncread(coverMapFilename, 'lon');
indexLatstartLC = find(latLandCoverMap<=maxLatArea,1, 'first');
indexLatendLC = find(latLandCoverMap<=minLatArea,1, 'first');
indexLonstartLC = find(lonLandCoverMap>=minLonArea,1, 'first');
indexLonendLC = find(lonLandCoverMap>=maxLonArea,1, 'first');
countLatLC=indexLatendLC-indexLatstartLC+1;
countLonLC=indexLonendLC-indexLonstartLC+1;
latLandCover_res=latLandCoverMap(indexLatstartLC:indexLatendLC);
lonLandCover_res=lonLandCoverMap(indexLonstartLC:indexLonendLC);
%define the coordinate matrixes of the land cover map
landCoverMap(:,:,1)=repmat(latLandCover_res,1,length(lonLandCover_res)); %lat matrix
landCoverMap(:,:,2)=repmat(lonLandCover_res',length(latLandCover_res),1); %lon matrix
%read the land cover map of the selected region
coverMapOrig = ncread(coverMapFilename,'lccs_class',[indexLonstartLC indexLatstartLC],[countLonLC countLatLC])';

%% simplify the land cover map to a limited number of classes
% 0 water bodies
% 1 broadleaved forest
% 2 needleleaved forest
% 3 cropland/grassland/shrubland
% 4 bare soils
% 5 flooded forest
% 6 flooded shurbland
% 7 permanent ice (not included in the map, ice is considered bare soil)

coverMapOrig(coverMapOrig==0)=4;%no_data=bare soils
coverMapOrig(coverMapOrig==10)=3;%cropland_rainfed=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==11)=3;%cropland_rainfed_herbaceous_cover=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==12)=3;%cropland_rainfed_tree_or_shrub_cover=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==20)=3;%cropland_irrigated=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==30)=3;%mosaic_cropland=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==40)=3;%mosaic_natural_vegetation=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==50)=1;%tree_broadleaved_evergreen_closed_to_open=broadleaved forest
coverMapOrig(coverMapOrig==60)=1;%tree_broadleaved_deciduous_closed_to_open=broadleaved forest
coverMapOrig(coverMapOrig==61)=1;% tree_broadleaved_deciduous_closed=broadleaved forest
coverMapOrig(coverMapOrig==62)=1;% tree_broadleaved_deciduous_open=broadleaved forest
coverMapOrig(coverMapOrig==70)=2;%tree_needleleaved_evergreen_closed_to_open=needleleaved forest
coverMapOrig(coverMapOrig==71)=2;%tree_needleleaved_evergreen_closed=needleleaved forest
coverMapOrig(coverMapOrig==72)=2;%tree_needleleaved_evergreen_open=needleleaved forest
coverMapOrig(coverMapOrig==80)=2;%tree_needleleaved_deciduous_closed_to_open=needleleaved forest
coverMapOrig(coverMapOrig==81)=2;%tree_needleleaved_deciduous_closed=needleleaved forest
coverMapOrig(coverMapOrig==82)=2;% tree_needleleaved_deciduous_open=needleleaved forest
coverMapOrig(coverMapOrig==90)=1;% tree_mixed=broadleaved forest
coverMapOrig(coverMapOrig==100)=1;%mosaic_tree_and_shrub=broadleaved forest
coverMapOrig(coverMapOrig==110)=3;%mosaic_herbaceous=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==120)=3;%shrubland=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==121)=3;%shrubland_evergreen=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==122)=3;%shrubland_deciduous=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==-126)=3;%grassland=cropland/grassland/shrubland
coverMapOrig(coverMapOrig==-116)=4;%lichens_and_mosses=bare soils
coverMapOrig(coverMapOrig==-106)=4;%sparse_vegetation=bare soils
coverMapOrig(coverMapOrig==-105)=4;%sparse_tree=bare soils
coverMapOrig(coverMapOrig==-104)=4;%sparse_shrub=bare soils
coverMapOrig(coverMapOrig==-103)=4;% sparse_herbaceous=bare soils
coverMapOrig(coverMapOrig==-96)=5;% tree_cover_flooded_fresh_or_brakish_water=flooded forest
coverMapOrig(coverMapOrig==-86)=5;%tree_cover_flooded_saline_water=flooded forest
coverMapOrig(coverMapOrig==-76)=6;%shrub_or_herbaceous_cover_flooded=flooded shurbland
coverMapOrig(coverMapOrig==-66)=4;%urban=bare soils
coverMapOrig(coverMapOrig==-56)=4;%bare_areas=bare soils
coverMapOrig(coverMapOrig==-55)=4;%bare_areas_consolidated=bare soils
coverMapOrig(coverMapOrig==-54)=4;%bare_areas_unconsolidated=bare soils
coverMapOrig(coverMapOrig==-46)=0;%water
coverMapOrig(coverMapOrig==-36)=4;%snow_and_ice=bare soils

%number of classes in the study area
coverMap_IDs=unique(coverMapOrig);
totNumOfCoverages=size(coverMap_IDs,1);
for l=1:totNumOfCoverages
    percentCoverageOfAclass(l)=length(coverMapOrig(coverMapOrig==coverMap_IDs(l)))/numel(coverMapOrig)*100;
end
maxPercentCoverageOfAclass=max(percentCoverageOfAclass);

% added by amir, change dominant class if it is waterbody to the next
% maximum (always dominant class is over land and not waterbody)
if percentCoverageOfAclass(coverMap_IDs == 0) == max(percentCoverageOfAclass)
    percentCoverageOfAclass_temp = percentCoverageOfAclass;
    percentCoverageOfAclass_temp(1) = 0;
    maxPercentCoverageOfAclass=max(percentCoverageOfAclass_temp);
    change_dominant_to_land_flag = 1;
end

dominantClass_ID=coverMap_IDs(percentCoverageOfAclass==maxPercentCoverageOfAclass);
for l=1:totNumOfCoverages
    if coverMap_IDs(l)~=0 && percentCoverageOfAclass(l)<0.1
        coverMapOrig(coverMapOrig==coverMap_IDs(l))=dominantClass_ID;
    end
    if coverMap_IDs(l)==0 && percentCoverageOfAclass(l)<0.01
        coverMapOrig(coverMapOrig==coverMap_IDs(l))=dominantClass_ID;
    end    
end

landCoverMap(:,:,3)=coverMapOrig;

if ~exist('change_dominant_to_land_flag', 'var')
    change_dominant_to_land_flag = 0;
end

end