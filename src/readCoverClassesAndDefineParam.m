%% ************************************************************************
% MODULE NAME:      readCoverClassesAndDefineParam.m
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
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [landCoverMap,coverMap_IDs,terrainParameters,vegetationParameters,landCoverMap_all]=readCoverClassesAndDefineParam(tag,coverMapFilename, SP_lla_onEllipsoid,DEMinfo,bioGeoInputs,presentLogInputsFilename)

if  tag.directsignalonly == 0
[landCoverMap,dominantClass_ID,change_dominant_to_land_flag]=readCoverMapCCI(coverMapFilename, SP_lla_onEllipsoid, DEMinfo.sizeInputMapWindow);

landCoverMap_all = landCoverMap;
%show land cover map in lat-lon
if tag.GUI ==  0
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

else
landCoverMap = 36 + (38-36) * rand(721, 721, 3);
landCoverMap(:,:,3) = 4;
dominantClass_ID = 0;
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


%list the found cover classes and ask the user to keep or change the
%parameters
if tag.GUI == 0
    if totNumOfCoverages>1
        %previous msg
        % inhomogSimul=questdlg([int2str(totNumOfCoverages) ' cover classes have been found in the simulation area.' ...                       
        %     'The class IDs are: ' int2str(coverMap_IDs')...
        %     'The class labels are: ' coverMap_label ...
        %     'Would you like to simulate all classes or only the dominant class?'], 'Coverage selection', 'All classes', 'Dominant class only', 'All classes');
        
        %% new msg, added by amir, to notify user about changing dominant
        %% class if the dominant class is waterbody
        if iscell(coverMap_label)
            labelsStr = strjoin(coverMap_label, sprintf('\n'));
        else
            labelsStr = coverMap_label;
        end

        idsStr = strtrim(int2str(coverMap_IDs'));

        % Notification depending on dominant class
        if change_dominant_to_land_flag == 1
            change_dominant_notif = sprintf('Dominant class was waterbody; changed to %d for land simulation.\n', dominantClass_ID);
        else
            change_dominant_notif = sprintf('Dominant class is: %d\n', dominantClass_ID);
        end

        msg = sprintf(['%d cover classes have been found in the simulation area.\n' ...
                       '%s' ...                          
                       'The class IDs are: %s\n' ...
                       'The class labels are:\n%s\n' ...
                       'Would you like to simulate all classes or only the dominant class?'], ...
                       totNumOfCoverages, change_dominant_notif, idsStr, labelsStr);

        inhomogSimul = questdlg(msg, ...
        'Coverage selection', ...
        'All classes', ... 
        'Dominant class only', ...
        'Change dominant class', ... 
        'All classes');


        switch inhomogSimul
            case 'Dominant class only'
                coverMap_label=coverMap_label(coverMap_IDs==dominantClass_ID);
                coverMap_IDs=dominantClass_ID;
                totNumOfCoverages=1;
                landCoverMap(:,:,3)=dominantClass_ID;
            %% Added by amir for changing dominant class to a specific class
            case 'Change dominant class'
                if iscell(coverMap_label)
                    labelsCell = coverMap_label(:);
                else
                    labelsCell = cellstr(coverMap_label);
                end
            
                ids = coverMap_IDs(:);
            
                n = min(numel(ids), numel(labelsCell));
                lines = cell(n,1);
                for k = 1:n
                    lines{k} = sprintf('%d: %s', ids(k), labelsCell{k});
                end
            
                listText = strjoin(lines, sprintf('\n'));
            
                prompt = sprintf("Select a new dominant class:\n\n%s\n\nEnter a class ID:", listText);
                answer = inputdlg(prompt, 'Change dominant class', [1 40]);
            
                if isempty(answer)
                    return;
                end
            
                New_selected_dominant_class = str2double(answer{1});
            
                if ~ismember(New_selected_dominant_class, coverMap_IDs)
                    errordlg('Invalid class ID.', 'Error');
                    return;
                end
            
                dominantClass_ID = New_selected_dominant_class;
                coverMap_label = coverMap_label(coverMap_IDs == dominantClass_ID);
                coverMap_IDs = dominantClass_ID;
                totNumOfCoverages = 1;
                landCoverMap(:,:,3) = dominantClass_ID;


        end
    end
else
    inhomogSimul='All classes'; % added by ansha
end


if tag.GUI == 0

    if totNumOfCoverages>1
        defaultCoverClass=questdlg([int2str(totNumOfCoverages) ' cover classes will be simulated.' ...                       
            'The class IDs are: ' int2str(coverMap_IDs')...
            'The class labels are: ' coverMap_label ...
            'Would you like to use the default inputs?'], 'Selection of terrain and vegetation inputs', 'yes', 'no', 'yes');
    else
       defaultCoverClass=questdlg(['One cover class will be simulated.' ...                       
            'The class IDs is: ' int2str(coverMap_IDs')...
            'The class label is: ' coverMap_label ...
            'Would you like to use the default inputs?'], 'Selection of terrain and vegetation inputs', 'yes', 'no', 'yes');  
    end
else
    defaultCoverClass='no'; % added by ansha
end


%% ---- default inputs ----
%flag for soil roughness autocorrelation function
%type of autocorrelation function: 0=gaussian 1D, 1=exponential 1D,
%2=gaussian 2D, 3=exponential 2D, 4=mixed-exponential 2D, 5=1.5-power
corr_funct = 3;

%flag for DBH distribution
%0=Log-normal, 1=Weibull
kindDBH_dist=0;

%% ---- open GUI for each class -----
switch defaultCoverClass
    case 'no'  
        
        for i=1:totNumOfCoverages   

            switch coverMap_IDs(i)
        
                case 0 %water bodies
                    if tag.GUI == 0
                        prompt={'Water status: Normal/Frozen',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Water temperature [ °C ]  :'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Normal','0.1', '20.0'};
                        defaultanswer=bioGeoInputs.water; % to check by ansha - 20230601 for no gui 
                       
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<3
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.water=classAnswer'; %values written in the previous inputs file
    
                        buttonSoilStatus=classAnswer{1};
                        switch buttonSoilStatus
                            case 'Normal'
                                mv_ini(i) = 80.0;   
                            case 'Frozen'
                                mv_ini(i) = 10.0;
                        end
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        surfaceTemperature(i)=str2double(classAnswer{3});
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        bulkDensity(i)=1.18;
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    else
                        defaultanswer=bioGeoInputs.water; % to check by ansha - 20230601 for no gui 
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.water=defaultanswer'; %values written in the previous inputs file
    
                        buttonSoilStatus=defaultanswer{1};
                        switch buttonSoilStatus
                            case 'Normal'
                                mv_ini(i) = 80.0;   
                            case 'Frozen'
                                mv_ini(i) = 10.0;
                        end
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        surfaceTemperature(i)=str2double(defaultanswer{3});
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        bulkDensity(i)=1.18;
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    end

                case 1 %broadleaved forest
                    if tag.GUI==0
                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ]  (set 80 for flooded case):',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:',...
                            'Total above ground biomass [t/ha]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Normal','1.0','5','30','1.18', '20.0', '5.','300.','25.','10.', '150.'};
                        defaultanswer=bioGeoInputs.BLForest; %values written in the previous inputs file   % to check by ansha - 20230601 for no gui
                       
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<11
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.BLForest=classAnswer'; % to check by ansha - 20230601 for no gui - skipped.
    
                        %load all values in the GUI
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6}); %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= str2num(classAnswer{11});   % double = above ground biomass        
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        if Bio(i) == 0.0 && LAI(i)>0.0
                            LAI(i) = 0.0;
                            CreateStruct.Interpreter = 'tex';
                            CreateStruct.WindowStyle = 'modal';
                            uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                                'LAI will be forced to zero.'},CreateStruct))
                        end
                else 
                    defaultanswer=bioGeoInputs.BLForest; %values written in the previous inputs file   % to check by ansha - 20230601 for no gui
                    % update the info to be written in the UserInputs.mat
                    bioGeoInputs.BLForest=defaultanswer'; % to check by ansha - 20230601 for no gui - skipped.

                    %load all values in the GUI
                    buttonSoilStatus=defaultanswer{1};
                    sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                    l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                    mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                    bulkDensity(i)=str2double(defaultanswer{5});
                    surfaceTemperature(i)=str2double(defaultanswer{6}); %the routine needs to have this input, but it is not used for the case of normal terrain
                    LAI(i)= str2num(defaultanswer{7});    % double = LAI
                    NTree(i)= str2num(defaultanswer{8});   % double = plant density
                    dbh_mean(i)= str2num(defaultanswer{9});   % double = mean dbh
                    DbhStD(i)= str2num(defaultanswer{10});   % double = std dbh
                    Bio(i)= str2num(defaultanswer{11});   % double = above ground biomass        
                    indexPlant(i) = 1 ;
                    kindPlant(i) = 1;
                    allomEq(i) = 0 ;
                    amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                    if Bio(i) == 0.0 && LAI(i)>0.0
                        LAI(i) = 0.0;
                        CreateStruct.Interpreter = 'tex';
                        CreateStruct.WindowStyle = 'modal';
                        uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                            'LAI will be forced to zero.'},CreateStruct))
                    end
            end
                case 2 %needleleaved forest
                    if tag.GUI ==0
                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ]  (set 80 for flooded case) :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:',...
                            'Total above ground biomass [t/ha]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Normal','1.0','5','15.0','1.18','20.0','5.','600.','15.','10.', '150.'};
                        defaultanswer=bioGeoInputs.NLForest; %values written in the previous inputs file
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<11
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.NLForest=classAnswer';
    
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6});
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= str2num(classAnswer{11});   % double = above ground biomass
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        if Bio(i) == 0.0 && LAI(i)>0.0
                            LAI(i) = 0.0;
                            CreateStruct.Interpreter = 'tex';
                            CreateStruct.WindowStyle = 'modal';
                            uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                                'LAI will be forced to zero.'},CreateStruct))
                        end
                    else
                        defaultanswer=bioGeoInputs.NLForest; %values written in the previous inputs file
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.NLForest=defaultanswer';
    
                        buttonSoilStatus=defaultanswer{1};
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(defaultanswer{5});
                        surfaceTemperature(i)=str2double(defaultanswer{6});
                        LAI(i)= str2num(defaultanswer{7});    % double = LAI
                        NTree(i)= str2num(defaultanswer{8});   % double = plant density
                        dbh_mean(i)= str2num(defaultanswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(defaultanswer{10});   % double = std dbh
                        Bio(i)= str2num(defaultanswer{11});   % double = above ground biomass
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; 
                        if Bio(i) == 0.0 && LAI(i)>0.0
                            LAI(i) = 0.0;
                            CreateStruct.Interpreter = 'tex';
                            CreateStruct.WindowStyle = 'modal';
                            uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                                'LAI will be forced to zero.'},CreateStruct))
                        end
                    end

                case 3 %cropland/grassland/shrubland, not-flooded/flooded
                    if tag.GUI==0
                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ]  (set 80 for flooded case):',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :',...
                            'LAI :',...
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Normal','1.0','5','15','1.18','20.0','2.','0.75', '10.', '4.'};
                        defaultanswer=bioGeoInputs.LowVeg; %values written in the previous inputs file
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<10
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.LowVeg=classAnswer';
    
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6}); %the routine needs to have this input, but it is not used for the case of normal soil and water
                        LAI(i) = str2double(classAnswer{7});    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{8});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{9});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{10});    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                    else
                        defaultanswer=bioGeoInputs.LowVeg; %values written in the previous inputs file
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.LowVeg=defaultanswer';
    
                        buttonSoilStatus=defaultanswer{1};
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(defaultanswer{5});
                        surfaceTemperature(i)=str2double(defaultanswer{6}); %the routine needs to have this input, but it is not used for the case of normal soil and water
                        LAI(i) = str2double(defaultanswer{7});    % double = leaf area index
                        amoipl(i)  = str2double(defaultanswer{8});    % double = vegetation water content
                        wleaf(i)   = str2double(defaultanswer{9});    % double = Leaf_length
                        bleaf(i)   = str2double(defaultanswer{10});    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                    end


                case 4 %bare areas
                    if tag.GUI ==0

                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ] (set 80 for flooded case) :',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Normal','1.0','5','15','1.18','20.'};
                        defaultanswer=bioGeoInputs.bare; %values written in the previous inputs file
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<6
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.bare=classAnswer';
    
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6});                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    else
                        defaultanswer=bioGeoInputs.bare; %values written in the previous inputs file
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.bare=defaultanswer';
    
                        buttonSoilStatus=defaultanswer{1};
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(defaultanswer{5});
                        surfaceTemperature(i)=str2double(defaultanswer{6});                    
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                    end

                case 5 %broadleaved forest flooded
                    if tag.GUI ==0
                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ]  (set 80 for flooded case):',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :',...
                            'LAI :',...
                            'Tree Density [#/ha] :',...
                            'Mean DBH [cm]:',...
                            'Standard deviation DBH [cm]:',...
                            'Total above ground biomass [t/ha]:'};
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Flooded','0.1','5','80','1.18', '20.0', '5.','300.','25.','10.', '150.'};
                        defaultanswer=bioGeoInputs.floodedForest; %values written in the previous inputs file
                       
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<11
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.floodedForest=classAnswer';
    
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6}); %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(classAnswer{7});    % double = LAI
                        NTree(i)= str2num(classAnswer{8});   % double = plant density
                        dbh_mean(i)= str2num(classAnswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(classAnswer{10});   % double = std dbh
                        Bio(i)= str2num(classAnswer{11});   % double = above ground biomass        
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        if Bio(i) == 0.0 && LAI(i)>0.0
                            LAI(i) = 0.0;
                            CreateStruct.Interpreter = 'tex';
                            CreateStruct.WindowStyle = 'modal';
                            uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                                'LAI will be forced to zero.'},CreateStruct))
                        end
                    else 
                       defaultanswer=bioGeoInputs.floodedForest; %values written in the previous inputs file
                       bioGeoInputs.floodedForest=defaultanswer';
                        buttonSoilStatus=defaultanswer{1};
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(defaultanswer{5});
                        surfaceTemperature(i)=str2double(defaultanswer{6}); %the routine needs to have this input, but it is not used for the case of normal terrain
                        LAI(i)= str2num(defaultanswer{7});    % double = LAI
                        NTree(i)= str2num(defaultanswer{8});   % double = plant density
                        dbh_mean(i)= str2num(defaultanswer{9});   % double = mean dbh
                        DbhStD(i)= str2num(defaultanswer{10});   % double = std dbh
                        Bio(i)= str2num(defaultanswer{11});   % double = above ground biomass        
                        indexPlant(i) = 1 ;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        if Bio(i) == 0.0 && LAI(i)>0.0
                            LAI(i) = 0.0;
                            CreateStruct.Interpreter = 'tex';
                            CreateStruct.WindowStyle = 'modal';
                            uiwait(msgbox({'Warning: Unfeasable condition in forest class: biomass is zero and LAI is larger than zero.';...
                                'LAI will be forced to zero.'},CreateStruct))
                        end
                    end


                case 6 %cropland/grassland/shrubland flooded
                    if tag.GUI ==0

                        prompt={'Surface status: Normal/Frozen/Flooded',...
                            'Soil roughness standard deviation [ cm ] :',...
                            'Soil roughness correlation length [ cm ] :',...
                            'Volumetric soil moisture [ % ]  (set 80 for flooded case):',...
                            'Dry soil bulk density [ g/cm^3 ]  :',...
                            'Surface temperature [ °C ]  :',...
                            'LAI :',...
                            'Plant gravimetric moisture content [g/g]:',...
                            'Leaf length [ cm ] :',...
                            'Leaf width [ cm ] :'};            
                        opts.Resize='on';
                        opts.WindowStyle='normal';
                        opts.Interpreter='tex';
                        name=['Class ' num2str(coverMap_IDs(i)) '  ' char(coverMap_label(i))];
                        numlines=[1,60]; %changed by ansha
                        %defaultanswer={'Flooded','0.1','5','80','1.18','20.0','2.','0.75', '10.', '4.'};
                        defaultanswer=bioGeoInputs.floodedLowVeg; %values written in the previous inputs file
            
                        classAnswer=inputdlg(prompt,name,numlines,defaultanswer,opts);
                        if size(classAnswer,1)<10
                            errordlg('Class inputs not found','Input error');
                            return
                        end
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.floodedLowVeg=classAnswer';
    
                        buttonSoilStatus=classAnswer{1};
                        sds_ini(i) = str2double(classAnswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(classAnswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(classAnswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(classAnswer{5});
                        surfaceTemperature(i)=str2double(classAnswer{6}); %the routine needs to have this input, but it is not used for the case of normal soil and water
                        LAI(i) = str2double(classAnswer{7});    % double = leaf area index
                        amoipl(i)  = str2double(classAnswer{8});    % double = vegetation water content
                        wleaf(i)   = str2double(classAnswer{9});    % double = Leaf_length
                        bleaf(i)   = str2double(classAnswer{10});    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                    else
                         defaultanswer=bioGeoInputs.floodedLowVeg; %values written in the previous inputs file
                        % update the info to be written in the UserInputs.mat
                        bioGeoInputs.floodedLowVeg=defaultanswer';
    
                        buttonSoilStatus=defaultanswer{1};
                        sds_ini(i) = str2double(defaultanswer{2});       % double = Rough_std.dev
                        l_ini(i) = str2double(defaultanswer{3}); % double = Rough.correl.L
                        mv_ini(i) = str2double(defaultanswer{4});   % double = Soil moisture
                        bulkDensity(i)=str2double(defaultanswer{5});
                        surfaceTemperature(i)=str2double(defaultanswer{6}); %the routine needs to have this input, but it is not used for the case of normal soil and water
                        LAI(i) = str2double(defaultanswer{7});    % double = leaf area index
                        amoipl(i)  = str2double(defaultanswer{8});    % double = vegetation water content
                        wleaf(i)   = str2double(defaultanswer{9});    % double = Leaf_length
                        bleaf(i)   = str2double(defaultanswer{10});    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                    end


            end
            switch buttonSoilStatus
                case 'Normal'
                    if coverMap_IDs(i)==0
                        terrainFlag(i)=2;
                    else
                        terrainFlag(i)=0;
                    end
                case 'Frozen'
                    terrainFlag(i)=1;
                case 'Flooded'
                    terrainFlag(i)=2;                            
            end
            terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
            vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                  Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];         
        end            
 %% ---- use default parameters ----
    case 'yes'
        for i=1:totNumOfCoverages
            switch coverMap_IDs(i)
                case 0
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
                        bioGeoInputs.water={'Normal','0.1', '20.0'};
                case 1
                        % 1 broadleaved forest
                        sds_ini(i) = 1.0;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 30.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=0;
                        LAI(i)= 5;    % double = LAI
                        NTree(i)= 300;   % double = plant density
                        dbh_mean(i)= 25;   % double = mean dbh
                        DbhStD(i)= 10;   % double = std dbh
                        Bio(i)= 150;   % double = above ground biomass
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        indexPlant(i) = 1;
                        bioGeoInputs.BLForest   = {'Normal','1.0','5.0','30','1.18', '20.0', '5.0','300.0','25.0','10.0', '150.0'};
                case 2
                        % 2 needleleaved forest
                        sds_ini(i) = 1.0;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 15.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=0;
                        LAI(i)= 5;    % double = LAI
                        NTree(i)= 600;   % double = plant density
                        dbh_mean(i)= 15;   % double = mean dbh
                        DbhStD(i)= 10;   % double = std dbh
                        Bio(i)= 150;   % double = above ground biomass
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        kindPlant(i) = 1;
                        allomEq(i) = 2 ;
                        indexPlant(i) = 1;
                        bioGeoInputs.NLForest   = {'Normal','1.0','5.0','15.0','1.18','20.0','5.0','600.','15.0','10.0', '150.0'};
                case 3
                        % 3 cropland/grassland/shrubland
                        sds_ini(i) = 1.0;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 15.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=0;
                        LAI(i) = 2;    % double = leaf area index
                        amoipl(i)  = 0.75;    % double = plant gravimetric moisture content
                        wleaf(i)   = 10.;    % double = Leaf_length
                        bleaf(i)   = 4.;    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0; allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        bioGeoInputs.LowVeg     = {'Normal','1.0','5.0','15.0','1.18','20.0','2.0','0.75', '10.0', '4.0'};
                case 4
                        % 4 bare soils
                        sds_ini(i) = 1.0;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 15.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=0;
                        indexPlant(i) = 0 ;
                        LAI(i)=0; amoipl(i)=0; wleaf(i)=0; bleaf(i)=0; kindPlant(i)=0;
                        Bio(i)=0; NTree(i)=0; dbh_mean(i)=0; DbhStD(i)=0; allomEq(i)=0;
                        bioGeoInputs.bare       = {'Normal','1.0','5.0','15.0','1.18','20.0'};
                case 5
                        % 5 flooded forest
                        sds_ini(i) = 0.1;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 80.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=2;
                        LAI(i)= 5;    % double = LAI
                        NTree(i)= 300;   % double = plant density
                        dbh_mean(i)= 25;   % double = mean dbh
                        DbhStD(i)= 10;   % double = std dbh
                        Bio(i)= 150;   % double = above ground biomass
                        amoipl(i)=0; wleaf(i)=0; bleaf(i)=0;
                        kindPlant(i) = 1;
                        allomEq(i) = 0 ;
                        indexPlant(i) = 1;
                        bioGeoInputs.floodedForest  = {'Flooded','0.1','5.0','80.0','1.18', '20.0', '5.0','300.','25.0','10.0', '150.0'};
                case 6
                        % 6 flooded shurbland                        
                        sds_ini(i) = 0.1;       % double = Rough_std.dev
                        l_ini(i) = 5.0; % double = Rough.correl.L
                        mv_ini(i) = 80.0;   % double = Soil moisture
                        bulkDensity(i)=1.18;
                        surfaceTemperature(i)=20.0; %the routine needs to have this input, but it is not used for the case of water
                        terrainFlag(i)=2;
                        LAI(i) = 2;    % double = leaf area index
                        amoipl(i)  = 0.75;    % double = plant gravimetric moisture content
                        wleaf(i)   = 10.;    % double = Leaf_length
                        bleaf(i)   = 4.;    % double = Leaf_width
                        Bio(i)=0;NTree(i)=0; dbh_mean(i)=0;DbhStD(i)=0;allomEq(i)=0;
                        kindPlant(i) = 0;
                        indexPlant(i) = 1;
                        bioGeoInputs.floodedLowVeg  = {'Flooded','0.1','5.0','80.0','1.18','20.0','2.0','0.75', '10.0', '4.0'};
            end           
                terrainParameters(i,:)=[terrainFlag(i) sds_ini(i) l_ini(i) corr_funct mv_ini(i) bulkDensity(i) surfaceTemperature(i)];
                vegetationParameters(i,:)=[indexPlant(i) kindPlant(i) LAI(i) amoipl(i) wleaf(i) bleaf(i) ...
                                            Bio(i) NTree(i) dbh_mean(i) DbhStD(i) kindDBH_dist allomEq(i)];
        end
end 
if tag.directsignalonly == 1;
else
if tag.GUI ==0  % added by ansha
  close(figLandCoverMap)
end
end
save('../conf/UserInputs.mat', 'bioGeoInputs', '-append');
save(presentLogInputsFilename, 'bioGeoInputs', '-append');
end