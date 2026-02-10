%% ************************************************************************
% MODULE NAME:      readCoverClassesAndDefineVariableParam.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente
% @copyright:   Tor Vergata University of Rome
% Original version:      2021 - 2022 by Laura Dente
% Main updates:          2022 by L. Dente
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module calls the routine to read the CCI land cover map and ask the
% user to provide the bio geo parameters for each cover class.
% This module is specific for the case of variable bio geo paremeters
% (those parameters varies in descrete steps or randomly).
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [landCoverMap,coverMap_IDs,terrainParameters,vegetationParameters,descreteVariableParam,randomVariableParameter,landCoverMap_all]=readCoverClassesAndDefineVariableParam(tag,coverMapFilename, SP_lla_onEllipsoid,DEMinfo,scenario, bioGeoInputsVariable, presentLogInputsFilename)


[landCoverMap,dominantClass_ID]=readCoverMapCCI(coverMapFilename, SP_lla_onEllipsoid, DEMinfo.sizeInputMapWindow);

landCoverMap_all = landCoverMap;

if tag.GUI == 0
%show land cover map in lat-lon
            mapOfColorCover=[0 0 1;0 1 0;1 0 1;0 1 1;1 0 0;1 1 0;0 0 0];
            %show land cover map in lat-lon
            %plot the SPs over the land cover
            figLandCoverMap = figure('Name', 'Land cover map', 'NumberTitle','off','OuterPosition', [50 50 700 510]);
            subplot('Position', [0.1 0.12 0.65 0.82]);
            mapshow(landCoverMap(:,:,2),landCoverMap(:,:,1),landCoverMap(:,:,3), 'DisplayType', 'texturemap')
%            axis([minLonArea maxLonArea minLatArea maxLatArea])
            colormap(mapOfColorCover)
            colorbar('Ticks',[0,1,2,3,4,5,6],...
                'TickLabels',{'0 water bodies','1 broadleaved forest','2 needleleaved forest','3 cropland/grassland/shrubland',...
                '4 bare areas','5 flooded forest','6 flooded shrubland'})
            caxis([0 6])
            title('Land cover map')
            xlabel('Longitude (deg)')
            ylabel('Latitude (deg)')
end
%% ----define soil and vegetation parameters for each land cover class

%number of classes in the study area
coverMap_IDs=unique(landCoverMap(:,:,3));
totNumOfCoverages=size(coverMap_IDs,1);
for l=1:totNumOfCoverages
    if coverMap_IDs(l)==0
        coverMap_label(l)={'water bodies'};
    elseif coverMap_IDs(l)==1
        coverMap_label(l)={'broadleaved forest'};
    elseif coverMap_IDs(l)==2
        coverMap_label(l)={'needleleaved forest'};
    elseif coverMap_IDs(l)==3
        coverMap_label(l)={'cropland/grassland/shrubland'};
    elseif coverMap_IDs(l)==4
        coverMap_label(l)={'bare areas'};
    elseif coverMap_IDs(l)==5
        coverMap_label(l)={'flooded forest'};
    elseif coverMap_IDs(l)==6
        coverMap_label(l)={'flooded shrubland'};
    end
        
end


%list the founded cover classes and ask the user to keep or change the
%parameters
if tag.GUI == 0
    if totNumOfCoverages>1
         inhomogSimul=questdlg([int2str(totNumOfCoverages) ' cover classes have been found in the simulation area.' ...                       
           'The class IDs are: ' int2str(coverMap_IDs')...
           'The class labels are: ' coverMap_label ...
           'Would you like to simulate all classes or only the dominant class?'], 'Coverage selection', 'All classes', 'Dominant class only', 'All classes');

         switch inhomogSimul
            case 'Dominant class only'
             coverMap_label=coverMap_label(coverMap_IDs==dominantClass_ID);
             coverMap_IDs=dominantClass_ID;
             totNumOfCoverages=1;
             landCoverMap(:,:,3)=dominantClass_ID;
     end
    end
 else
    inhomogSimul='All classes'; % added by ansha
end

%flag for soil roughness autocorrelation function
%type of autocorrelation function: 0=gaussian 1D, 1=exponential 1D,
%2=gaussian 2D, 3=exponential 2D, 4=mixed-exponential 2D, 5=1.5-power
corr_funct = 3;

%flag for DBH distribution
%0=Log-normal, 1=Weibull
kindDBH_dist=0;

switch scenario 
    case 1 %soil moisture scenario
       if tag.GUI == 0
            prompt={'Minimum volumetric soil moisture [ % ] :',...
                    'Maximum volumetric soil moisture [ % ] :',...
                    'Step volumetric soil moisture [ % ] :',...
                    'Number of random variations of the other variables for each soil moisture value (max 9) [ # ] :'};
            opts.Resize='on';
            opts.WindowStyle='normal';
            opts.Interpreter='tex';
            name='Inputs for the descrete variability of soil moisture ';
            numlines=[1,70]; %changed by ansha
            %defaultanswer={'1.0','41.0','10','3'};
            defaultanswer=bioGeoInputsVariable.sm.inputsVar; %get the inputs from the previous inputs file
    
            smAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
            if size(smAnswer,1)<4
                errordlg('Incorrect inputs for the soil moisture variability','Input error');
                return
            end
            bioGeoInputsVariable.sm.inputsVar=smAnswer'; %update the UserInputs file and log file
            mv_min = str2double(smAnswer{1}); % double = Soil moisture min
            mv_max = str2double(smAnswer{2}); % double = Soil moisture max
            mv_step = str2double(smAnswer{3}); % double = Soil moisture increment
            numRuns = str2double(smAnswer{4}); % double = numbers of runs to repeat for each soil moisture value
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[mv_min mv_step mv_max numRuns];
       else
            smAnswer = bioGeoInputsVariable.sm.inputsVar;
            bioGeoInputsVariable.sm.inputsVar=smAnswer'; %update the UserInputs file and log file
            mv_min = str2double(smAnswer{1}); % double = Soil moisture min
            mv_max = str2double(smAnswer{2}); % double = Soil moisture max
            mv_step = str2double(smAnswer{3}); % double = Soil moisture increment
            numRuns = str2double(smAnswer{4}); % double = numbers of runs to repeat for each soil moisture value
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[mv_min mv_step mv_max numRuns];
       end

        
        for i=1:totNumOfCoverages   

            switch coverMap_IDs(i)
        
                case 0 %water bodies
                    %water bodies
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = 5.0; % double = Rough.correl.L
                    mv_ini(i) = 80.0;   % double = Soil moisture
                    bulkDensity(i)=1.18;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                    terrainFlag(i)=2;
                    indexPlant(i) = 0 ;
                    LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                    Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    rough_min(i)=sds_ini(i); rough_max(i)=rough_min(i);
                    bio_min(i)=0; bio_max(i)=0; LAI_min(i)=0; LAI_max(i)=0;

                case 1 %broadleaved forest
                    if tag.GUI == 0 
                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','25.','275.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.BLForest{1},bioGeoInputsVariable.sm.BLForest{2},...
                                    bioGeoInputsVariable.sm.BLForest{3},bioGeoInputsVariable.sm.BLForest{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.sm.BLForest{5},bioGeoInputsVariable.sm.BLForest{6},...
                                    bioGeoInputsVariable.sm.BLForest{7},bioGeoInputsVariable.sm.BLForest{8},...
                                    bioGeoInputsVariable.sm.BLForest{9},bioGeoInputsVariable.sm.BLForest{10}};
                       
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<12
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.BLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        bio_min(i)=str2num(classAnswer{4});
                        bio_max(i)=str2num(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        LAI_min(i)=0; LAI_max(i)=0;
                    else                       
                      %  classAnswer=bioGeoInputsVariable.sm.BLForest; % by
                      %  Mauro to align with case of GUI where element 1
                      %  and 6 of prompt are empty
                        classAnswer=['---', bioGeoInputsVariable.sm.BLForest(1:4) '....' bioGeoInputsVariable.sm.BLForest(5:10) ] ; 
            
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.BLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        bio_min(i)=str2num(classAnswer{4});
                        bio_max(i)=str2num(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        LAI_min(i)=0; LAI_max(i)=0;
                   end
                
                case 2 %needleleaved forest
                    if tag.GUI == 0
                    prompt={'---------------------------------------------------------------------',...
                        'Random soil roughness standard deviation: minimum value [ cm ] :',...
                        'Random soil roughness standard deviation: maximum value [ cm ] :',...
                        'Random above ground biomass: minimum value [t/ha]:',...
                        'Random above ground biomass: maximum value [t/ha]:',...
                        '---------------------------------------------------------------------',...
                        'Soil roughness correlation length [ cm ] :',...
                        'Dry soil bulk density [ g/cm^3 ]  :',...
                        'LAI :',...
                        'Tree Density [#/ha] :',...
                        'Mean DBH [cm]:',...
                        'Standard deviation DBH [cm]:'};
                    opts.Resize='on';
                    opts.WindowStyle='normal';
                    opts.Interpreter='tex';
                    name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                    numlines=[1,70]; %changed by ansha
                    %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0', '25.','275',...
                    %                '           FIXED PARAMETERS','5','1.18','5.','600.','15.','10.'};
                    %get the inputs from the previous inputs file
                    defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.NLForest{1},bioGeoInputsVariable.sm.NLForest{2},...
                                bioGeoInputsVariable.sm.NLForest{3},bioGeoInputsVariable.sm.NLForest{4},...
                                '           FIXED PARAMETERS',bioGeoInputsVariable.sm.NLForest{5},bioGeoInputsVariable.sm.NLForest{6},...
                                bioGeoInputsVariable.sm.NLForest{7},bioGeoInputsVariable.sm.NLForest{8},...
                                bioGeoInputsVariable.sm.NLForest{9},bioGeoInputsVariable.sm.NLForest{10}};

                    classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                    if size(classAnswer,1)<12
                        errordlg('Class inputs not found','Input error');
                        return
                    end
                    %update the UserInputs file and log file
                    bioGeoInputsVariable.sm.NLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                        classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};                    
                    rough_min(i) = str2double(classAnswer{2});
                    rough_max(i) = str2double(classAnswer{3});
                    bio_min(i)=str2num(classAnswer{4});
                    bio_max(i)=str2num(classAnswer{5});
                    sds_ini(i) = 0.0;       % double = Rough_std.dev
                    l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                    mv_ini(i) = 0.0;   % double = Soil moisture
                    bulkDensity(i)=str2double(classAnswer{8});
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                    LAI(i)= str2num(classAnswer{9});    % double = LAI
                    NTree(i)= str2num(classAnswer{10});   % double = plant density
                    dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                    DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                    Bio(i)= 0.0;   % double = above ground biomass                    
                    terrainFlag(i)=0;    
                    indexPlant(i) = 1 ;
                    kindPlant(i) = 1;
                    allomEq(i) = 2 ;
                    amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                    LAI_min(i)=0; LAI_max(i)=0;        
                    else
                        % changed by mauro to mimic the case with GUI where
                        % promt ha size 12
                     classAnswer=['....', bioGeoInputsVariable.sm.NLForest(1:4) '....' bioGeoInputsVariable.sm.NLForest(5:10) ];
         
                    %update the UserInputs file and log file
                    bioGeoInputsVariable.sm.NLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                        classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};                    
                    rough_min(i) = str2double(classAnswer{2});
                    rough_max(i) = str2double(classAnswer{3});
                    bio_min(i)=str2num(classAnswer{4});
                    bio_max(i)=str2num(classAnswer{5});
                    sds_ini(i) = 0.0;       % double = Rough_std.dev
                    l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                    mv_ini(i) = 0.0;   % double = Soil moisture
                    bulkDensity(i)=str2double(classAnswer{8});
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                    LAI(i)= str2num(classAnswer{9});    % double = LAI
                    NTree(i)= str2num(classAnswer{10});   % double = plant density
                    dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                    DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                    Bio(i)= 0.0;   % double = above ground biomass                    
                    terrainFlag(i)=0;    
                    indexPlant(i) = 1 ;
                    kindPlant(i) = 1;
                    allomEq(i) = 2 ;
                    amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                    LAI_min(i)=0; LAI_max(i)=0; 
                    end

                case 3 %cropland/grassland/shrubland
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.LowVeg{1},bioGeoInputsVariable.sm.LowVeg{2},...
                                    bioGeoInputsVariable.sm.LowVeg{3},bioGeoInputsVariable.sm.LowVeg{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.sm.LowVeg{5},bioGeoInputsVariable.sm.LowVeg{6},...
                                    bioGeoInputsVariable.sm.LowVeg{7},bioGeoInputsVariable.sm.LowVeg{8},...
                                    bioGeoInputsVariable.sm.LowVeg{9}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<11
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.LowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        LAI_min(i) = str2double(classAnswer{4});
                        LAI_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});                  
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0; bio_min(i)=0;bio_max(i)=0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                        % modified by mauri
                        classAnswer= ['....' bioGeoInputsVariable.sm.LowVeg(1:4) '....' bioGeoInputsVariable.sm.LowVeg(5:9) ];
        
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.LowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        LAI_min(i) = str2double(classAnswer{4});
                        LAI_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});                  
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0; bio_min(i)=0;bio_max(i)=0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    end

                case 4 %bare areas
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0',...
                        %                '           FIXED PARAMETERS','5','1.18'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.bare{1},bioGeoInputsVariable.sm.bare{2},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.sm.bare{3},bioGeoInputsVariable.sm.bare{4}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<6
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.bare={classAnswer{2},classAnswer{3},classAnswer{5},classAnswer{6}};
                        
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{6});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        bio_min(i)=0;bio_max(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.sm.bare;
    
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.bare={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4}};
                        
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{4});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        bio_min(i)=0;bio_max(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    end

                case 5 % flooded broadleaved forest
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','25.','275.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.floodedForest{1},bioGeoInputsVariable.sm.floodedForest{2},...
                                    bioGeoInputsVariable.sm.floodedForest{3},bioGeoInputsVariable.sm.floodedForest{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.sm.floodedForest{5},bioGeoInputsVariable.sm.floodedForest{6},...
                                    bioGeoInputsVariable.sm.floodedForest{7},bioGeoInputsVariable.sm.floodedForest{8},...
                                    bioGeoInputsVariable.sm.floodedForest{9},bioGeoInputsVariable.sm.floodedForest{10}};
                       
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<12
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        bio_min(i)=str2num(classAnswer{4});
                        bio_max(i)=str2num(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        LAI_min(i)=0; LAI_max(i)=0;
                    else
                        classAnswer=['....', bioGeoInputsVariable.sm.floodedForest(1:4) '....' bioGeoInputsVariable.sm.floodedForest(5:10) ];
    
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        bio_min(i)=str2num(classAnswer{4});
                        bio_max(i)=str2num(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        LAI_min(i)=0; LAI_max(i)=0;
                    end

                 case 6 %cropland/grassland/shrubland, flooded
                     if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.sm.floodedLowVeg{1},bioGeoInputsVariable.sm.floodedLowVeg{2},...
                                    bioGeoInputsVariable.sm.floodedLowVeg{3},bioGeoInputsVariable.sm.floodedLowVeg{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.sm.floodedLowVeg{5},bioGeoInputsVariable.sm.floodedLowVeg{6},...
                                    bioGeoInputsVariable.sm.floodedLowVeg{7},bioGeoInputsVariable.sm.floodedLowVeg{8},...
                                    bioGeoInputsVariable.sm.floodedLowVeg{9}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<11
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.floodedLowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        LAI_min(i) = str2double(classAnswer{4});
                        LAI_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});                  
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0; bio_min(i)=0;bio_max(i)=0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;     
                     else
                        classAnswer=['....' bioGeoInputsVariable.sm.floodedLowVeg(1:4) '....' bioGeoInputsVariable.sm.floodedLowVeg(5:9) ] ;

                        %update the UserInputs file and log file
                        bioGeoInputsVariable.sm.floodedLowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        LAI_min(i) = str2double(classAnswer{4});
                        LAI_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{8});                  
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0; bio_min(i)=0;bio_max(i)=0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;  
                     end


            end
            terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
            vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                  Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];
            randomVariableParameter(i,:)=[rough_min(i) rough_max(i) bio_min(i) bio_max(i) LAI_min(i) LAI_max(i)];
            
        end
    case 2 %biomass scenario
        if tag.GUI ==0 
            prompt={'Minimum above ground biomass [ t/ha ] :',...
                    'Maximum above ground biomass [ t/ha ] :',...
                    'Step above ground biomass [ t/ha ] :',...
                    'Number of random variations of the other variables for each biomass value (max 9) [ # ] :'};
            opts.Resize='on';
            opts.WindowStyle='normal';
            opts.Interpreter='tex';
            name='Inputs for the descrete variability of above ground biomass ';
            numlines=[1,70]; %changed by ansha
            %defaultanswer={'25.0','275.0','50','3'};
            defaultanswer=bioGeoInputsVariable.bio.inputsVar; %get the inputs from the previous inputs file
    
            bioAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
            if size(bioAnswer,1)<4
                errordlg('Incorrect inputs for the biomass variability','Input error');
                return
            end
            bioGeoInputsVariable.bio.inputsVar=bioAnswer';%update the UserInputs file and log file
    
            min_bio = str2double(bioAnswer{1}); % double = Soil moisture min
            max_bio = str2double(bioAnswer{2}); % double = Soil moisture max
            bio_step = str2double(bioAnswer{3}); % double = Soil moisture increment
            numRuns = str2double(bioAnswer{4}); % double = numbers of runs to repeat for each biomass value
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[min_bio bio_step max_bio numRuns];
        else
            bioAnswer=bioGeoInputsVariable.bio.inputsVar; %get the inputs from the previous inputs file
    
            bioGeoInputsVariable.bio.inputsVar=bioAnswer';%update the UserInputs file and log file
    
            min_bio = str2double(bioAnswer{1}); % double = Soil moisture min
            max_bio = str2double(bioAnswer{2}); % double = Soil moisture max
            bio_step = str2double(bioAnswer{3}); % double = Soil moisture increment
            numRuns = str2double(bioAnswer{4}); % double = numbers of runs to repeat for each biomass value
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[min_bio bio_step max_bio numRuns];
        end

        
        for i=1:totNumOfCoverages   

            switch coverMap_IDs(i)
        
                case 0 %water bodies
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = 5.0; % double = Rough.correl.L
                    mv_ini(i) = 80.0;   % double = Soil moisture
                    bulkDensity(i)=1.18;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                    terrainFlag(i)=2;
                    indexPlant(i) = 0 ;
                    LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                    Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    rough_min(i)=sds_ini(i); rough_max(i)=sds_ini(i);moist_min(i)=mv_ini(i); moist_max(i)=mv_ini(i);
                    LAI_min(i)=0; LAI_max(i)=0;

                case 1 %broadleaved forest
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.BLForest{1},bioGeoInputsVariable.bio.BLForest{2},...
                                    bioGeoInputsVariable.bio.BLForest{3},bioGeoInputsVariable.bio.BLForest{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.BLForest{5},bioGeoInputsVariable.bio.BLForest{6},...
                                    bioGeoInputsVariable.bio.BLForest{7},bioGeoInputsVariable.bio.BLForest{8},...
                                    bioGeoInputsVariable.bio.BLForest{9},bioGeoInputsVariable.bio.BLForest{10}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<12
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.BLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};                    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else 
                         classAnswer=bioGeoInputsVariable.bio.BLForest ;
  
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.BLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},...
                                            classAnswer{5},classAnswer{6},classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10}};                    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{6});
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    end

                
                case 2 %needleleaved forest
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','600.','15.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.NLForest{1},bioGeoInputsVariable.bio.NLForest{2},...
                                    bioGeoInputsVariable.bio.NLForest{3},bioGeoInputsVariable.bio.NLForest{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.NLForest{5},bioGeoInputsVariable.bio.NLForest{6},...
                                    bioGeoInputsVariable.bio.NLForest{7},bioGeoInputsVariable.bio.NLForest{8},...
                                    bioGeoInputsVariable.bio.NLForest{9},bioGeoInputsVariable.bio.NLForest{10}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<12
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.NLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;        
                    else
                        classAnswer= bioGeoInputsVariable.bio.NLForest;
   
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.NLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},...
                                            classAnswer{5},classAnswer{6},classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{6});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0; 
                    end

                           
                case 3 %cropland/grassland/shrubland
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.LowVeg{1},bioGeoInputsVariable.bio.LowVeg{2},...
                                    bioGeoInputsVariable.bio.LowVeg{3},bioGeoInputsVariable.bio.LowVeg{4},...
                                    bioGeoInputsVariable.bio.LowVeg{5},bioGeoInputsVariable.bio.LowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.LowVeg{7},bioGeoInputsVariable.bio.LowVeg{8},...
                                    bioGeoInputsVariable.bio.LowVeg{9},bioGeoInputsVariable.bio.LowVeg{10},...
                                    bioGeoInputsVariable.bio.LowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.LowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.bio.LowVeg;

                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.LowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        LAI_min(i) = str2double(classAnswer{5});
                        LAI_max(i) = str2double(classAnswer{6}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    end


                case 4 %bare areas
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.bare{1},bioGeoInputsVariable.bio.bare{2},...
                                    bioGeoInputsVariable.bio.bare{3},bioGeoInputsVariable.bio.bare{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.bare{5},bioGeoInputsVariable.bio.bare{6}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<8
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.bare={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8}};
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    else
                         classAnswer=bioGeoInputsVariable.bio.bare;
         
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.bare={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},...
                                            classAnswer{5},classAnswer{6}};
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{6});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    end


                case 5 %broadleaved forest flooded
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.floodedForest{1},bioGeoInputsVariable.bio.floodedForest{2},...
                                    bioGeoInputsVariable.bio.floodedForest{3},bioGeoInputsVariable.bio.floodedForest{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.floodedForest{5},bioGeoInputsVariable.bio.floodedForest{6},...
                                    bioGeoInputsVariable.bio.floodedForest{7},bioGeoInputsVariable.bio.floodedForest{8},...
                                    bioGeoInputsVariable.bio.floodedForest{9},bioGeoInputsVariable.bio.floodedForest{10}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<12
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};                    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.bio.floodedForest;
        
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.floodedForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},...
                                            classAnswer{5},classAnswer{6},classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10}};                    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{6});
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    end


                case 6 %cropland/grassland/shrubland, flooded
                    if tag.GUI == 0
                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.bio.floodedLowVeg{1},bioGeoInputsVariable.bio.floodedLowVeg{2},...
                                    bioGeoInputsVariable.bio.floodedLowVeg{3},bioGeoInputsVariable.bio.floodedLowVeg{4},...
                                    bioGeoInputsVariable.bio.floodedLowVeg{5},bioGeoInputsVariable.bio.floodedLowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.bio.floodedLowVeg{7},bioGeoInputsVariable.bio.floodedLowVeg{8},...
                                    bioGeoInputsVariable.bio.floodedLowVeg{9},bioGeoInputsVariable.bio.floodedLowVeg{10},...
                                    bioGeoInputsVariable.bio.floodedLowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.floodedLowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.bio.floodedLowVeg;

                        %update the UserInputs file and log file
                        bioGeoInputsVariable.bio.floodedLowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        LAI_min(i) = str2double(classAnswer{5});
                        LAI_max(i) = str2double(classAnswer{6}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;

                    end


            end
            terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
            vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                  Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];
            randomVariableParameter(i,:)=[rough_min(i) rough_max(i) moist_min(i) moist_max(i) LAI_min(i) LAI_max(i)];
            
        end

    case 3 %freeze-thaw scenario
        if tag.GUI == 0

            prompt={'Number of runs with random bio/geo parameters (max 9) [ # ] :'};
            opts.Resize='on';
            opts.WindowStyle='normal';
            opts.Interpreter='tex';
            name='Number of simulation repetition for the freeze-thaw scenario ';
            numlines=[1,70];
            %defaultanswer={'3'};
            defaultanswer=bioGeoInputsVariable.ft.inputsVar; %get the inputs from the previous inputs file
    
            numberRunAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
            if size(numberRunAnswer,1)<1
                errordlg('Incorrect inputs for the number of Runs','Input error');
                return
            end
            bioGeoInputsVariable.ft.inputsVar=numberRunAnswer;%update the UserInputs file and log file
    
            numRuns = str2double(numberRunAnswer{1}); % double = numbers of runs to repeat (in total will be done x2 runs, first unfrozen then frozen)
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[0 1 1 numRuns];
        else
             defaultanswer=bioGeoInputsVariable.ft.inputsVar; %get the inputs from the previous inputs file
    
            numberRunAnswer=defaultanswer;
            if size(numberRunAnswer,1)<1
                errordlg('Incorrect inputs for the number of Runs','Input error');
                return
            end
            bioGeoInputsVariable.ft.inputsVar=numberRunAnswer;%update the UserInputs file and log file
    
            numRuns = str2double(numberRunAnswer{1}); % double = numbers of runs to repeat (in total will be done x2 runs, first unfrozen then frozen)
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[0 1 1 numRuns];
        end


        for i=1:totNumOfCoverages   

            switch coverMap_IDs(i)
        
                case 0 %water bodies
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = 5.0; % double = Rough.correl.L
                    mv_ini(i) = 80.0;   % double = Soil moisture
                    bulkDensity(i)=1.18;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                    terrainFlag(i)=2;
                    indexPlant(i) = 0 ;
                    LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                    bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    rough_min(i)=sds_ini(i); rough_max(i)=sds_ini(i);moist_min(i)=mv_ini(i); moist_max(i)=mv_ini(i);
                    LAI_min(i)=0; LAI_max(i)=0;

                case 1 %broadleaved forest
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275.',...
                        %               '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'}; 
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.BLForest{1},bioGeoInputsVariable.ft.BLForest{2},...
                                    bioGeoInputsVariable.ft.BLForest{3},bioGeoInputsVariable.ft.BLForest{4},...
                                    bioGeoInputsVariable.ft.BLForest{5},bioGeoInputsVariable.ft.BLForest{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.BLForest{7},bioGeoInputsVariable.ft.BLForest{8},...
                                    bioGeoInputsVariable.ft.BLForest{9},bioGeoInputsVariable.ft.BLForest{10},...
                                    bioGeoInputsVariable.ft.BLForest{11},bioGeoInputsVariable.ft.BLForest{12}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<14
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.BLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        bio_min(i)=str2num(classAnswer{6});
                        bio_max(i)=str2num(classAnswer{7});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});
                        LAI(i)= str2num(classAnswer{11});    % double = LAI
                        NTree(i)= str2num(classAnswer{12});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{14});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else
                       classAnswer=bioGeoInputsVariable.ft.BLForest;
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.BLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        bio_min(i)=str2num(classAnswer{5});
                        bio_max(i)=str2num(classAnswer{6});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    end

                    
                case 2 %needleleaved forest
                    if tag.GUI == 0 

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275.',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','600.','15.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.NLForest{1},bioGeoInputsVariable.ft.NLForest{2},...
                                    bioGeoInputsVariable.ft.NLForest{3},bioGeoInputsVariable.ft.NLForest{4},...
                                    bioGeoInputsVariable.ft.NLForest{5},bioGeoInputsVariable.ft.NLForest{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.NLForest{7},bioGeoInputsVariable.ft.NLForest{8},...
                                    bioGeoInputsVariable.ft.NLForest{9},bioGeoInputsVariable.ft.NLForest{10},...
                                    bioGeoInputsVariable.ft.NLForest{11},bioGeoInputsVariable.ft.NLForest{12}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<14
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.NLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        bio_min(i)=str2num(classAnswer{6});
                        bio_max(i)=str2num(classAnswer{7}); 
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{11});    % double = LAI
                        NTree(i)= str2num(classAnswer{12});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{14});   % double = std dbh                  
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;     
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.ft.NLForest;
               
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.NLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        bio_min(i)=str2num(classAnswer{5});
                        bio_max(i)=str2num(classAnswer{6}); 
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh                  
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;     
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    end
                
                case 3 %cropland/grassland/shrubland, not-flooded
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.LowVeg{1},bioGeoInputsVariable.ft.LowVeg{2},...
                                    bioGeoInputsVariable.ft.LowVeg{3},bioGeoInputsVariable.ft.LowVeg{4},...
                                    bioGeoInputsVariable.ft.LowVeg{5},bioGeoInputsVariable.ft.LowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.LowVeg{7},bioGeoInputsVariable.ft.LowVeg{8},...
                                    bioGeoInputsVariable.ft.LowVeg{9},bioGeoInputsVariable.ft.LowVeg{10},...
                                    bioGeoInputsVariable.ft.LowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.LowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.ft.LowVeg;
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.LowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        LAI_min(i) = str2double(classAnswer{5});
                        LAI_max(i) = str2double(classAnswer{6}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    end


                case 4 %bare areas
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.bare{1},bioGeoInputsVariable.ft.bare{2},...
                                    bioGeoInputsVariable.ft.bare{3},bioGeoInputsVariable.ft.bare{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.bare{5},bioGeoInputsVariable.ft.bare{6}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<8
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.bare={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8}}; 
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.ft.bare;

                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.bare={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},...
                                            classAnswer{5},classAnswer{6}}; 
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{5}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{6});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;

                    end


                case 5 %flooded forest (needleleaved)
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275.',...
                        %               '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'}; 
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.floodedForest{1},bioGeoInputsVariable.ft.floodedForest{2},...
                                    bioGeoInputsVariable.ft.floodedForest{3},bioGeoInputsVariable.ft.floodedForest{4},...
                                    bioGeoInputsVariable.ft.floodedForest{5},bioGeoInputsVariable.ft.floodedForest{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.floodedForest{7},bioGeoInputsVariable.ft.floodedForest{8},...
                                    bioGeoInputsVariable.ft.floodedForest{9},bioGeoInputsVariable.ft.floodedForest{10},...
                                    bioGeoInputsVariable.ft.floodedForest{11},bioGeoInputsVariable.ft.floodedForest{12}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<14
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        bio_min(i)=str2num(classAnswer{6});
                        bio_max(i)=str2num(classAnswer{7});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});
                        LAI(i)= str2num(classAnswer{11});    % double = LAI
                        NTree(i)= str2num(classAnswer{12});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{14});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ; %(needledleaved)
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else
                      classAnswer=bioGeoInputsVariable.ft.floodedForest;
                       
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.floodedForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        bio_min(i)=str2num(classAnswer{5});
                        bio_max(i)=str2num(classAnswer{6});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ; %(needledleaved)
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    end


                case 6 %cropland/grassland/shrubland, flooded
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75', '10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.ft.floodedLowVeg{1},bioGeoInputsVariable.ft.floodedLowVeg{2},...
                                    bioGeoInputsVariable.ft.floodedLowVeg{3},bioGeoInputsVariable.ft.floodedLowVeg{4},...
                                    bioGeoInputsVariable.ft.floodedLowVeg{5},bioGeoInputsVariable.ft.floodedLowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.ft.floodedLowVeg{7},bioGeoInputsVariable.ft.floodedLowVeg{8},...
                                    bioGeoInputsVariable.ft.floodedLowVeg{9},bioGeoInputsVariable.ft.floodedLowVeg{10},...
                                    bioGeoInputsVariable.ft.floodedLowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.floodedLowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                       classAnswer=bioGeoInputsVariable.ft.floodedLowVeg;
          
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.ft.floodedLowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        LAI_min(i) = str2double(classAnswer{5});
                        LAI_max(i) = str2double(classAnswer{6}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    end


            end
            terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
            vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                  Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];
            randomVariableParameter(i,:)=[rough_min(i) rough_max(i) moist_min(i) moist_max(i) LAI_min(i) LAI_max(i) bio_min(i) bio_max(i)];
            
        end        
    case 4 %wetland
        if tag.GUI == 0 

            prompt={'Number of runs with random bio/geo parameters [ # ] :'};
            opts.Resize='on';
            opts.WindowStyle='normal';
            opts.Interpreter='tex';
            name='Number of simulation repetition for the wetland scenario (max 9) ';
            numlines=[1,70];
            %defaultanswer={'3'};
            defaultanswer=bioGeoInputsVariable.wet.inputsVar; %get the inputs from the previous inputs file        
    
            numberRunAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
            if size(numberRunAnswer,1)<1
                errordlg('Incorrect inputs for the number of Runs','Input error');
                return
            end
            bioGeoInputsVariable.wet.inputsVar=numberRunAnswer;%update the UserInputs file and log file
    
            numRuns = str2double(numberRunAnswer{1}); % double = numbers of runs to repeat (in total will be done x2 runs, first normal then flooded)
                if numRuns>9
                    errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                    numRuns = 9;
                end
                descreteVariableParam=[0 2 2 numRuns];
        else
        %get the inputs from the previous inputs file        
    
            numberRunAnswer=bioGeoInputsVariable.wet.inputsVar;

            bioGeoInputsVariable.wet.inputsVar=numberRunAnswer;%update the UserInputs file and log file
    
            numRuns = str2double(numberRunAnswer{1}); % double = numbers of runs to repeat (in total will be done x2 runs, first normal then flooded)
            if numRuns>9
                errordlg('The number of simulation repetitions is higher than 9. It will be reset to 9.','Input error');
                numRuns = 9;
            end
            descreteVariableParam=[0 2 2 numRuns];    
        end


        for i=1:totNumOfCoverages   

            switch coverMap_IDs(i)
        
                case 0 %water bodies
                    %water bodies
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = 5.0; % double = Rough.correl.L
                    mv_ini(i) = 80.0;   % double = Soil moisture
                    bulkDensity(i)=1.18;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                    terrainFlag(i)=2;
                    indexPlant(i) = 0 ;
                    LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                    bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    rough_min(i)=sds_ini(i); rough_max(i)=sds_ini(i);moist_min(i)=mv_ini(i); moist_max(i)=mv_ini(i);
                    LAI_min(i)=0; LAI_max(i)=0;

                case 1 %broadleaved forest
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275',...
                        %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.BLForest{1},bioGeoInputsVariable.wet.BLForest{2},...
                                    bioGeoInputsVariable.wet.BLForest{3},bioGeoInputsVariable.wet.BLForest{4},...
                                    bioGeoInputsVariable.wet.BLForest{5},bioGeoInputsVariable.wet.BLForest{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.wet.BLForest{7},bioGeoInputsVariable.wet.BLForest{8},...
                                    bioGeoInputsVariable.wet.BLForest{9},bioGeoInputsVariable.wet.BLForest{10},...
                                    bioGeoInputsVariable.wet.BLForest{11},bioGeoInputsVariable.wet.BLForest{12}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<14
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.BLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        bio_min(i)=str2num(classAnswer{6});
                        bio_max(i)=str2num(classAnswer{7});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});
                        LAI(i)= str2num(classAnswer{11});    % double = LAI
                        NTree(i)= str2num(classAnswer{12});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{14});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.wet.BLForest;
               
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.BLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        bio_min(i)=str2num(classAnswer{5});
                        bio_max(i)=str2num(classAnswer{6});
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh
                        Bio(i)= 0.0;   % double = above ground biomass
                        terrainFlag(i)=0;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;                        
                    end

                
                case 2 %needleleaved forest
                    if tag.GUI == 0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random above ground biomass: minimum value [t/ha]:',...
                            'Random above ground biomass: maximum value [t/ha]:',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275.',...
                        %               '           FIXED PARAMETERS','5','1.18','5.','600.','15.','10.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.NLForest{1},bioGeoInputsVariable.wet.NLForest{2},...
                                    bioGeoInputsVariable.wet.NLForest{3},bioGeoInputsVariable.wet.NLForest{4},...
                                    bioGeoInputsVariable.wet.NLForest{5},bioGeoInputsVariable.wet.NLForest{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.wet.NLForest{7},bioGeoInputsVariable.wet.NLForest{8},...
                                    bioGeoInputsVariable.wet.NLForest{9},bioGeoInputsVariable.wet.NLForest{10},...
                                    bioGeoInputsVariable.wet.NLForest{11},bioGeoInputsVariable.wet.NLForest{12}};
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<14
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.NLForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        bio_min(i)=str2num(classAnswer{6});
                        bio_max(i)=str2num(classAnswer{7}); 
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{11});    % double = LAI
                        NTree(i)= str2num(classAnswer{12});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{14});   % double = std dbh                  
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;     
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;   
                    else
                        classAnswer= bioGeoInputsVariable.wet.NLForest;
    
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.NLForest={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        bio_min(i)=str2num(classAnswer{5});
                        bio_max(i)=str2num(classAnswer{6}); 
                        sds_ini(i) = 0.0;       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{9});    % double = LAI
                        NTree(i)= str2num(classAnswer{10});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{11});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{12});   % double = std dbh                  
                        Bio(i)= 0.0;   % double = above ground biomass                   
                        terrainFlag(i)=0;     
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;  
                    end

            
                case 3 %cropland/grassland/shrubland, not-flooded
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75','10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.LowVeg{1},bioGeoInputsVariable.wet.LowVeg{2},...
                                    bioGeoInputsVariable.wet.LowVeg{3},bioGeoInputsVariable.wet.LowVeg{4},...
                                    bioGeoInputsVariable.wet.LowVeg{5},bioGeoInputsVariable.wet.LowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.wet.LowVeg{7},bioGeoInputsVariable.wet.LowVeg{8},...
                                    bioGeoInputsVariable.wet.LowVeg{9},bioGeoInputsVariable.wet.LowVeg{10},...
                                    bioGeoInputsVariable.wet.LowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.LowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{10});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    else
                        classAnswer=bioGeoInputsVariable.wet.LowVeg;
        
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.LowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                            classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};
    
                        rough_min(i) = str2double(classAnswer{1});
                        rough_max(i) = str2double(classAnswer{2});
                        moist_min(i) = str2double(classAnswer{3});
                        moist_max(i) = str2double(classAnswer{4});
                        LAI_min(i) = str2double(classAnswer{5});
                        LAI_max(i) = str2double(classAnswer{6}); 
                        sds_ini(i) = 0.0; 
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});                   
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=0;
                    end


                case 4 %bare areas
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation: minimum value [ cm ] :',...
                            'Random soil roughness standard deviation: maximum value [ cm ] :',...
                            'Random volumetric soil moisture: minimum value [ % ] :',...
                            'Random volumetric soil moisture: maximum value [ % ] :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.',...
                        %                '           FIXED PARAMETERS','5','1.18'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.bare{1},bioGeoInputsVariable.wet.bare{2},...
                                    bioGeoInputsVariable.wet.bare{3},bioGeoInputsVariable.wet.bare{4},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.wet.bare{5},bioGeoInputsVariable.wet.bare{6}};
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<8
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.bare={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8}}; 
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    else
                        classAnswer= ['....' bioGeoInputsVariable.wet.bare(1:4) '....' bioGeoInputsVariable.wet.bare(5:6) ];
                  
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.bare={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},...
                                            classAnswer{7},classAnswer{8}}; 
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        sds_ini(i) = 0.0;                    
                        l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                        mv_ini(i) = 0.0;
                        bulkDensity(i)=str2double(classAnswer{8});
                        surfaceTemperature(i)=20.0;                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        LAI_min(i)=0; LAI_max(i)=0;
                        terrainFlag(i)=0;
                    end


                case 5 %broadleaved forest flooded
                    if tag.GUI ==0
                        prompt={'---------------------------------------------------------------------',...
                        'Random soil roughness standard deviation : minimum value [ cm ] :',...
                        'Random soil roughness standard deviation : maximum value [ cm ] :',...
                        'Random volumetric soil moisture (when not flooded): minimum value [ % ] :',...
                        'Random volumetric soil moisture (when not flooded): maximum value [ % ] :',...
                        'Random above ground biomass: minimum value [t/ha]:',...
                        'Random above ground biomass: maximum value [t/ha]:',...
                        '---------------------------------------------------------------------',...
                        'Soil roughness correlation length [ cm ] :',...
                        'Dry soil bulk density [ g/cm^3 ]  :',...
                        'LAI :',...
                        'Tree Density [#/ha] :',...
                        'Mean DBH [cm]:',...
                        'Standard deviation DBH [cm]:'};
                    opts.Resize='on';
                    opts.WindowStyle='normal';
                    opts.Interpreter='tex';
                    name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                    numlines=[1,70]; %changed by ansha
                    %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','25.','275.',...
                    %                '           FIXED PARAMETERS','5','1.18','5.','300.','25.','10.'};
                    %get the inputs from the previous inputs file
                    defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.floodedForest{1},bioGeoInputsVariable.wet.floodedForest{2},...
                                bioGeoInputsVariable.wet.floodedForest{3},bioGeoInputsVariable.wet.floodedForest{4},...
                                bioGeoInputsVariable.wet.floodedForest{5},bioGeoInputsVariable.wet.floodedForest{6},...
                                '           FIXED PARAMETERS',bioGeoInputsVariable.wet.floodedForest{7},bioGeoInputsVariable.wet.floodedForest{8},...
                                bioGeoInputsVariable.wet.floodedForest{9},bioGeoInputsVariable.wet.floodedForest{10},...
                                bioGeoInputsVariable.wet.floodedForest{11},bioGeoInputsVariable.wet.floodedForest{12}};

                    classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                    if size(classAnswer,1)<14
                        errordlg('Class inputs not found','Input error');
                        return
                    end
                    %update the UserInputs file and log file
                    bioGeoInputsVariable.wet.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                        classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};

                    rough_min(i) = str2double(classAnswer{2});
                    rough_max(i) = str2double(classAnswer{3});
                    moist_min(i) = str2double(classAnswer{4});
                    moist_max(i) = str2double(classAnswer{5});
                    bio_min(i)=str2num(classAnswer{6});
                    bio_max(i)=str2num(classAnswer{7});
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                    mv_ini(i) = 99.0;
                    bulkDensity(i)=str2double(classAnswer{10});
                    LAI(i)= str2num(classAnswer{11});    % double = LAI
                    NTree(i)= str2num(classAnswer{12});   % double = plant density
                    dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                    DbhStD(i)= str2num(classAnswer{14});   % double = std dbh
                    Bio(i)= 0.0;   % double = above ground biomass
                    terrainFlag(i)=2;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                    indexPlant(i) = 1 ;
                    kindPlant(i) = 1;
                    allomEq(i) = 0 ;
                    amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                else
                    classAnswer=['....' bioGeoInputsVariable.wet.floodedForest(1:6) '....' bioGeoInputsVariable.wet.floodedForest(7:12) ];
         
                    %update the UserInputs file and log file
                    bioGeoInputsVariable.wet.floodedForest={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                        classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13},classAnswer{14}};

                    rough_min(i) = str2double(classAnswer{2});
                    rough_max(i) = str2double(classAnswer{3});
                    moist_min(i) = str2double(classAnswer{4});
                    moist_max(i) = str2double(classAnswer{5});
                    bio_min(i)=str2num(classAnswer{6});
                    bio_max(i)=str2num(classAnswer{7});
                    sds_ini(i) = 0.1;       % double = Rough_std.dev
                    l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                    mv_ini(i) = 99.0;
                    bulkDensity(i)=str2double(classAnswer{10});
                    LAI(i)= str2num(classAnswer{11});    % double = LAI
                    NTree(i)= str2num(classAnswer{12});   % double = plant density
                    dbh_mean(i)= str2num(classAnswer{13});   % double = mean dbh
                    DbhStD(i)= str2num(classAnswer{14});   % double = std dbh
                    Bio(i)= 0.0;   % double = above ground biomass
                    terrainFlag(i)=2;
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of normal terrain            
                    indexPlant(i) = 1 ;
                    kindPlant(i) = 1;
                    allomEq(i) = 0 ;
                    amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; LAI_min(i)=0; LAI_max(i)=0;
                end

                case 6 %cropland/grassland/shrubland flooded
                    if tag.GUI ==0

                        prompt={'---------------------------------------------------------------------',...
                            'Random soil roughness standard deviation (when not flooded): minimum value [ cm ] :',...
                            'Random soil roughness standard deviation (when not flooded): maximum value [ cm ] :',...
                            'Random volumetric soil moisture (when not flooded): minimum value [ % ] :',...
                            'Random volumetric soil moisture (when not flooded): maximum value [ % ] :',...
                            'Random LAI: minimum value :',...
                            'Random LAI: maximum value :',...
                            '---------------------------------------------------------------------',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...            
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,70]; %changed by ansha
                        %defaultanswer={'           RANDOMLY VARIABLE PARAMETERS','0.3','3.0','1.0','41.','0.','6.',...
                        %                '           FIXED PARAMETERS','5','1.18','0.75','10.', '4.'};
                        %get the inputs from the previous inputs file
                        defaultanswer={'           RANDOMLY VARIABLE PARAMETERS',bioGeoInputsVariable.wet.floodedLowVeg{1},bioGeoInputsVariable.wet.floodedLowVeg{2},...
                                    bioGeoInputsVariable.wet.floodedLowVeg{3},bioGeoInputsVariable.wet.floodedLowVeg{4},...
                                    bioGeoInputsVariable.wet.floodedLowVeg{5},bioGeoInputsVariable.wet.floodedLowVeg{6},...
                                    '           FIXED PARAMETERS',bioGeoInputsVariable.wet.floodedLowVeg{7},bioGeoInputsVariable.wet.floodedLowVeg{8},...
                                    bioGeoInputsVariable.wet.floodedLowVeg{9},bioGeoInputsVariable.wet.floodedLowVeg{10},...
                                    bioGeoInputsVariable.wet.floodedLowVeg{11}};
    
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<13
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        %update the UserInputs file and log file
                        bioGeoInputsVariable.wet.floodedLowVeg={classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},classAnswer{7},...
                                            classAnswer{9},classAnswer{10},classAnswer{11},classAnswer{12},classAnswer{13}};
    
                        rough_min(i) = str2double(classAnswer{2});
                        rough_max(i) = str2double(classAnswer{3});
                        moist_min(i) = str2double(classAnswer{4});
                        moist_max(i) = str2double(classAnswer{5});
                        LAI_min(i) = str2double(classAnswer{6});
                        LAI_max(i) = str2double(classAnswer{7});
                        sds_ini(i) = 0.1; 
                        l_ini(i) = str2double(classAnswer{9}); % double = Rough.correl.L
                        mv_ini(i) = 99.0;
                        bulkDensity(i)=str2double(classAnswer{10});                    
                        LAI(i) = 0.0;    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{11});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{12});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{13});    % double = Leaf_width
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        terrainFlag(i)=2;
                else

                    classAnswer= bioGeoInputsVariable.wet.floodedLowVeg;
       
                    %update the UserInputs file and log file
                    bioGeoInputsVariable.wet.floodedLowVeg={classAnswer{1},classAnswer{2},classAnswer{3},classAnswer{4},classAnswer{5},classAnswer{6},...
                                        classAnswer{7},classAnswer{8},classAnswer{9},classAnswer{10},classAnswer{11}};

                    rough_min(i) = str2double(classAnswer{1});
                    rough_max(i) = str2double(classAnswer{2});
                    moist_min(i) = str2double(classAnswer{3});
                    moist_max(i) = str2double(classAnswer{4});
                    LAI_min(i) = str2double(classAnswer{5});
                    LAI_max(i) = str2double(classAnswer{6});
                    sds_ini(i) = 0.1; 
                    l_ini(i) = str2double(classAnswer{7}); % double = Rough.correl.L
                    mv_ini(i) = 99.0;
                    bulkDensity(i)=str2double(classAnswer{8});                    
                    LAI(i) = 0.0;    % double = leaf area index
                    amoipl(i)  = str2double(classAnswer{9});    % double = vegetation water content
                    wleaf(i)   = str2double(classAnswer{10});    % double = Leaf_length
                    bleaf(i)   = str2double(classAnswer{11});    % double = Leaf_width
                    surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                    bio_min(i)=0; bio_max(i)=0; Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                    kindPlant(i) = 0;
                    indexPlant(i) = 1;
                    terrainFlag(i)=2;
               end

            end
            terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
            vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                  Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];
            randomVariableParameter(i,:)=[rough_min(i) rough_max(i) moist_min(i) moist_max(i) LAI_min(i) LAI_max(i) bio_min(i) bio_max(i)];
            
        end
end
if tag.GUI == 0    
    close(figLandCoverMap)
end
save('../conf/UserInputs.mat', 'bioGeoInputsVariable', '-append');
save(presentLogInputsFilename, 'bioGeoInputsVariable', '-append');
end