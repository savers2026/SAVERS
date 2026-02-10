%% ************************************************************************
% MODULE NAME:      computeGridCellParameters_dc.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, D. Comite, N. Pierdicca, L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      2010 by Roberto Giusto
% Main updates:          2018 by L. Dente (spaceborne)
%                        2021-2022 by D. Comite (HSAVERS, optimization, attitude)
%                        2022 by L. Dente (ionospheric delay)
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module computes all the parameters needed for the integration, for each facet
% of the simulation area: simulation grid, delay (including the ionospheric delay)
% doppler shift, local angles, antenna gain after pearcing.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [Xplot,Yplot,regGridSim_coverMap,sigmaParam,DDMparam,gammaParam] = computeGridCell_dc(TRAp, dataDEM_enuRot, coverMap_enuRot,...
                                                                                               geoSYSp, nameOutputFile, tag, DEMinfo, ionosphericDelay,coverMap_enuRot_all)

if tag.exe == 0, dbstop if error, end

% waitbarSimGrid = waitbar(0,'COMPUTATION OF SIMULATION GRID AND LOCAL ANGLES');

% pixels distances vs. GPS_Tx
sigmaParam.R0 = sqrt(TRAp.Txlocal'*TRAp.Txlocal);    

%% --------- ISODOPPLER and isodelay curves calculation -------------

% [m/usec] Light Speed per usec
speedlight = 2.99792458e2;

%vector difference Tx e Rx
RTx = TRAp.Txlocal - TRAp.Rxlocal;

%distance in [m] Tx - Rx (Euclidean norm or vector magnitude)
RTx_mod    = norm(RTx); 
Rx_mod     = norm(TRAp.Rxlocal);
Tx_mod     = norm(TRAp.Txlocal);

gammaParam.distanceRatio = (Rx_mod + Tx_mod)^2/RTx_mod^2;
gammaParam.squaredDistanceSum = (Rx_mod + Tx_mod)^2;
                        
% Reference delay in the rotated ENU origin
% the extra path introduced by the ionosphere along Tx-SP-Rx is included
refdelay = Rx_mod + Tx_mod + ionosphericDelay;
DDMparam.absoluteDelayAtSP_microsec = refdelay./speedlight;              % [microseconds: usec]
DDMparam.delayAtSPReferencedToDirect_microsec = (refdelay - RTx_mod)./speedlight;              % [microseconds: usec]

%% ------ definition of the simulation grid -------- 

%define the simulation grid in X (set the SP in the coordinate zero)
XplotPos = 0.0:geoSYSp.dx:max(geoSYSp.Xmax_Range);
XplotNeg = -(0.0+geoSYSp.dx:geoSYSp.dx:abs(min(geoSYSp.Xmin_Range)));
%local grid in X (array)
Xplot    = [fliplr(XplotNeg) XplotPos];
%number of grid cells in X
Nx_surf  = length(Xplot); 

%define the simulation grid in Y (set the SP in the coordinate zero)
YplotPos = 0.0:geoSYSp.dy:max(geoSYSp.Yabs_Range);
YplotNeg = -YplotPos(2:length(YplotPos));
Yplot    = [fliplr(YplotNeg) YplotPos];
%local grid in Y (array)
Yplot    = flipud(Yplot');
%number of grid cells in Y
Ny_surf  = length(Yplot);

% % Separation of simulation grid area for L1 and L5 (adeed by Amir)
% if geoSYSp.indexTxSignal == 1 || geoSYSp.indexTxSignal == 3
%     for indexTxFrequency = 1:length(TRAp.TxFrequency)
%         %define the simulation grid in X (set the SP in the coordinate zero)
%         XplotPos = 0.0:geoSYSp.dx:geoSYSp.Xmax_Range(indexTxFrequency);
%         XplotNeg = -(0.0+geoSYSp.dx:geoSYSp.dx:abs(geoSYSp.Xmin_Range(indexTxFrequency)));
%         %local grid in X (array)
%         Xplot{indexTxFrequency}    = [fliplr(XplotNeg) XplotPos];
%         %number of grid cells in X
%         Nx_surf(indexTxFrequency)  = length(Xplot{indexTxFrequency}); 
% 
%         %define the simulation grid in Y (set the SP in the coordinate zero)
%         YplotPos = 0.0:geoSYSp.dy:geoSYSp.Yabs_Range(indexTxFrequency);
%         YplotNeg = -YplotPos(2:length(YplotPos));
%         Yplot{indexTxFrequency}    = [fliplr(YplotNeg) YplotPos];
%         %local grid in Y (array)
%         Yplot{indexTxFrequency}    = flipud(Yplot{indexTxFrequency}');
%         %number of grid cells in Y
%         Ny_surf(indexTxFrequency)  = length(Yplot{indexTxFrequency});
%     end
% elseif geoSYSp.indexTxSignal == 0 || geoSYSp.indexTxSignal == 2
%     %define the simulation grid in X (set the SP in the coordinate zero)
%     XplotPos = 0.0:geoSYSp.dx:geoSYSp.Xmax_Range;
%     XplotNeg = -(0.0+geoSYSp.dx:geoSYSp.dx:abs(geoSYSp.Xmin_Range));
%     %local grid in X (array)
%     Xplot{1}    = [fliplr(XplotNeg) XplotPos];
%     %number of grid cells in X
%     Nx_surf(1)  = length(Xplot{1}); 
% 
%     %define the simulation grid in Y (set the SP in the coordinate zero)
%     YplotPos = 0.0:geoSYSp.dy:geoSYSp.Yabs_Range;
%     YplotNeg = -YplotPos(2:length(YplotPos));
%     Yplot{1}    = [fliplr(YplotNeg) YplotPos];
%     %local grid in Y (array)
%     Yplot{1}    = flipud(Yplot{1}');
%     %number of grid cells in Y
%     Ny_surf(1)  = length(Yplot{1});
% end
% 
% % end Separation of simulation grid area for L1 and L5 (adeed by Amir)

% regGridSim_x = cell(1, numel(Xplot));
% regGridSim_y = cell(1, numel(Yplot));
% 
% for i = 1:numel(TRAp.TxFrequency)
%     [regGridSim_x{i}, regGridSim_y{i}] = meshgrid(double(Xplot{i}), double(Yplot{i}));
% end

[regGridSim_x, regGridSim_y] = meshgrid(Xplot,Yplot);

%% --------- registration of DEM and land cover map over the simulation grid ------------------

%%%%%Note: the following griddata are time consuming

if DEMinfo.flatsurface == 0
    %find z component corresponding to the points of the
    %regular grid applying a linear interpolation of DEM
    regGridSim_z = griddata(dataDEM_enuRot(:,:,1), dataDEM_enuRot(:,:,2), dataDEM_enuRot(:,:,3), regGridSim_x, regGridSim_y,'linear');
else
    %set to zero the ENU z coordinates for a test without topography and
    %without surface curvature
    regGridSim_z = zeros(Nx_surf,Ny_surf);
end



% for i = 1:numel(TRAp.TxFrequency)
%     if DEMinfo.flatsurface == 0
%         %find z component corresponding to the points of the
%         %regular grid applying a linear interpolation of DEM
%         regGridSim_z{i} = griddata(dataDEM_enuRot(:,:,1), dataDEM_enuRot(:,:,2), dataDEM_enuRot(:,:,3), regGridSim_x{i}, regGridSim_y{i},'linear');
%     else
%         %set to zero the ENU z coordinates for a test without topography and
%         %without surface curvature
%         regGridSim_z{i} = zeros(Nx_surf(i),Ny_surf(i));
%     end
% end


%interpolate the land cover map over the simulation grid, using the nearest
%neighbour method
if length(unique(coverMap_enuRot(:,:,3)))>1
    %interpolate the land cover map over the simulation grid, using the nearest
    %neighbour method
    regGridSim_coverMap = griddata(coverMap_enuRot(:,:,1),coverMap_enuRot(:,:,2),coverMap_enuRot(:,:,3), regGridSim_x, regGridSim_y,'nearest');
else
    regGridSim_coverMap = zeros(Ny_surf,Nx_surf);
    regGridSim_coverMap(:,:) = coverMap_enuRot(1,1,3); % just this in original
end

regGridSim_coverMap_all = griddata(coverMap_enuRot_all(:,:,1),coverMap_enuRot_all(:,:,2),coverMap_enuRot_all(:,:,3), regGridSim_x, regGridSim_y,'nearest');


%% AMIR
regGridSim_z(regGridSim_coverMap_all == 0) = 0; % set the surface height to 0 if the surface is waterbody
% remove anomalies in height
globalAnomalyMask = abs((regGridSim_z - mean(regGridSim_z(:),'omitnan')) ./ std(regGridSim_z(:),'omitnan')) > 1;
filteredDEM = regGridSim_z; 
filteredDEM(globalAnomalyMask) = 0;
localDeviation = abs(filteredDEM - imfilter(filteredDEM, fspecial('average',[2 2]), 'replicate'));
filteredDEM(globalAnomalyMask | localDeviation > mean(localDeviation(:),'omitnan') + 2*std(localDeviation(:),'omitnan')) = 0;
regGridSim_z = filteredDEM;
%%


% regGridSim_coverMap = cell(1, numel(TRAp.TxFrequency));
% for i = 1:numel(TRAp.TxFrequency)
%     if length(unique(coverMap_enuRot(:,:,3)))>1
%         %interpolate the land cover map over the simulation grid, using the nearest
%         %neighbour method
%         regGridSim_coverMap{i} = griddata(coverMap_enuRot(:,:,1),coverMap_enuRot(:,:,2),coverMap_enuRot(:,:,3), regGridSim_x{i}, regGridSim_y{i},'nearest');
%     else
%         regGridSim_coverMap{i} = zeros(Ny_surf(i),Nx_surf(i));
%         regGridSim_coverMap{i}(:,:) = coverMap_enuRot(1,1,3);
%     end
% end

%Show and save DEM interpolated over the simulation grid
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
if tag.plotDEM == 1
    figDEMregGrid = figure('Name', 'DEM on simulation grid', 'NumberTitle','off');
    imagesc(Xplot,Yplot,regGridSim_z)
    colorbar, colormap(jet)
    title('DEM of the simulation area')
    axis xy
    xlabel('X (m)'), ylabel('Y (m)')
    saveas(figDEMregGrid, strjoin([nameOutputFile '_DEMonSimGrid.fig'],''));
    close(figDEMregGrid)
end
end

if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
if tag.plotCoverMap == 1
    %show cover map interpolated on the simulation regular grid  
    mapOfColorCover = [0 0 1;0 1 0;1 0 1;0 1 1;1 0 0;1 1 0;0 0 0];
    %show land cover map on the regular grid
    figLandCoverRegGrid=figure('Name', 'Land cover map on the simulation grid', 'NumberTitle','off','OuterPosition', [50 50 700 510]);
    subplot('Position', [0.1 0.12 0.65 0.82]);
    mapshow(regGridSim_x,regGridSim_y,regGridSim_coverMap, 'DisplayType', 'texturemap')
    colormap(mapOfColorCover)
    colorbar('Ticks',[0,1,2,3,4,5,6],...
                'TickLabels',{'0 water bodies','1 broadleaved forest','2 needleleaved forest','3 cropland/grassland/shrubland',...
                '4 bare areas','5 flooded forest','6 flooded shrubland'})
    caxis([0 6])
    title('Land cover map on the simulation grid')
    xlabel('X (m)')
    ylabel('Y (m)')
    saveas(figLandCoverRegGrid, strjoin([nameOutputFile '_LandCoverOnSimGrid.fig'],''));
%     if geoSYSp.Mode ~= "Single Point" , close(figLandCoverRegGrid), end ; 
end
end


% for i = 1:numel(TRAp.TxFrequency)
% 
%     %Show and save DEM interpolated over the simulation grid
%     if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
%     if tag.plotDEM == 1
%         figDEMregGrid = figure('Name', 'DEM on simulation grid', 'NumberTitle','off');
%         imagesc(Xplot{i},Yplot{i},regGridSim_z{i})
%         colorbar, colormap(jet)
%         title('DEM of the simulation area')
%         axis xy
%         xlabel('X (m)'), ylabel('Y (m)')
%         if numel(TRAp.TxFrequency) == 1
%             saveas(figDEMregGrid, strjoin([nameOutputFile '_DEMonSimGrid_L1.fig'],''));
%         elseif numel(TRAp.TxFrequency) == 2
%             if i == 1
%                 saveas(figDEMregGrid, strjoin([nameOutputFile '_DEMonSimGrid_L1.fig'],''));
%             elseif i == 2
%                 saveas(figDEMregGrid, strjoin([nameOutputFile '_DEMonSimGrid_L5.fig'],''));
%             end
%         end
% %         close(figDEMregGrid)
%     end
%     end
% 
%     if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1
%     if tag.plotCoverMap == 1
%         %show cover map interpolated on the simulation regular grid  
%         mapOfColorCover = [0 0 1;0 1 0;1 0 1;0 1 1;1 0 0;1 1 0;0 0 0];
%         %show land cover map on the regular grid
%         figLandCoverRegGrid=figure('Name', 'Land cover map on the simulation grid', 'NumberTitle','off','OuterPosition', [50 50 700 510]);
%         subplot('Position', [0.1 0.12 0.65 0.82]);
%         mapshow(regGridSim_x{i},regGridSim_y{i},regGridSim_coverMap{i}, 'DisplayType', 'texturemap')
%         colormap(mapOfColorCover)
%         colorbar('Ticks',[0,1,2,3,4,5,6],...
%                     'TickLabels',{'0 water bodies','1 broadleaved forest','2 needleleaved forest','3 cropland/grassland/shrubland',...
%                     '4 bare areas','5 flooded forest','6 flooded shrubland'})
%         caxis([0 6])
%         title('Land cover map on the simulation grid')
%         xlabel('X (m)')
%         ylabel('Y (m)')
%         if numel(TRAp.TxFrequency) == 1
%             saveas(figLandCoverRegGrid, strjoin([nameOutputFile '_LandCoverOnSimGrid_L1.fig'],''));
%         elseif numel(TRAp.TxFrequency) == 2
%             if i == 1
%                 saveas(figLandCoverRegGrid, strjoin([nameOutputFile '_LandCoverOnSimGrid_L1.fig'],''));
%             elseif i == 2
%                 saveas(figLandCoverRegGrid, strjoin([nameOutputFile '_LandCoverOnSimGrid_L5.fig'],''));
%             end
%         end
%     %     if geoSYSp.Mode ~= "Single Point" , close(figLandCoverRegGrid), end ; 
%     end
%     end
% 
% end


clear dataDEM_enuRot
clear coverMap_enuRot

%% -------- computation of slope and aspect at each grid cell --------
%compute slope and aspect of the interpolated DEM
[gradientX, gradientY] = gradient(regGridSim_z, geoSYSp.dx, -geoSYSp.dy);
[aspect, mag]          = cart2pol(gradientX, gradientY);

slope = atan(mag);

%set to NaN the aspect of the points with zero gradient
%this step is necessary to recognize these in points in the following steps
aspect(gradientY == 0 & gradientX == 0) = NaN;

aspect = wrapTo2Pi(-aspect-pi/2);

%transform aspect and slope in degree
[aspect, slope] = fromRadians('degrees', aspect, slope);

clear mag
clear gradientX
clear gradientY

%% ------- computation of dA and the factor of the integral equation ------------

%compute the effective surface areas of each grid cell, taking into account 
%of surface slope
[~, regGridSim_cellareas] = demarea(regGridSim_z, geoSYSp.dx, 'true');

DDMparam.regGridSim_cellareas = regGridSim_cellareas;

%factor used to compute the signal power of the down-looking antenna,
%without calibration (the values are more reasonable without Ti^2
for indexTxFrequency=1:length(TRAp.waveL)
    DDMparam.factor(:,:,indexTxFrequency) = TRAp.waveL(indexTxFrequency)^2*db2pow(TRAp.TxSignalPower(indexTxFrequency))*...
                                        TRAp.TxGainAtSP(indexTxFrequency)*regGridSim_cellareas./(4*pi)^3;
end

clear regGridSim_cellareas

%%   calculations for each Point (xp, yp) in the Area of Interest
%%  main loop computing for each facet of the grid (xp,yp,zp) the scattering angles,
%%  the delay and doppler, gain of Rx and Tx

DDMparam.fdoppler = ones(Ny_surf,Nx_surf,length(TRAp.waveL));     % matrix of doppler shifts vs. Reference
DDMparam.fdelay   = zeros(Ny_surf,Nx_surf,length(TRAp.waveL));    % matrix of range delay [usec] as required

sigmaParam.incAngle_flr        = zeros(Ny_surf,Nx_surf);
sigmaParam.incAzimuthAngle_flr = zeros(Ny_surf,Nx_surf);

sigmaParam.TH_s         = zeros(Ny_surf,Nx_surf);  %matrix scattering angle in Elevation, angle between scattering direction and normal to xy rot ENU plane
sigmaParam.PH_s         = zeros(Ny_surf,Nx_surf);  %matrix scattering angle in Azimuth
sigmaParam.TH_s_flr     = zeros(Ny_surf,Nx_surf);
sigmaParam.PH_s_flr     = zeros(Ny_surf,Nx_surf);
sigmaParam.deltaPhi_flr = zeros(Ny_surf,Nx_surf);

sigmaParam.zVect_flr    = zeros(Ny_surf,Nx_surf,3);
sigmaParam.r2mat        = ones(Ny_surf,Nx_surf);     % matrix of the Rx-to-pixels distances
DDMparam.PATH_factor    = ones(Ny_surf,Nx_surf);     % matrix of space attenuation (path lenght)

sigmaParam.prR_real     = zeros(Ny_surf,Nx_surf,2);
sigmaParam.prR_im       = zeros(Ny_surf,Nx_surf,2);
sigmaParam.prL_real     = zeros(Ny_surf,Nx_surf,2);
sigmaParam.prL_im       = zeros(Ny_surf,Nx_surf,2);

%Defining a new variable to project the Rx Down pattern on the area of interest
DDMparam.Gain_onsurf_lin = zeros(Ny_surf, Nx_surf,length(TRAp.waveL));
        
%%%% matrices of the coordinates
xp = regGridSim_x;
yp = regGridSim_y;
zp = regGridSim_z;

%% ------- DOWN LOOKING ANTENNA GAIN -------%
%---------------- SGR2ECEF -----------------%
        
%%% pSGR coordinate del punto di integrazione in SGR 
%%% SGR: SAVERS Global Reference, i.e., rotated enu  
pSGR(:,:,1) = xp;
pSGR(:,:,2) = yp;
pSGR(:,:,3) = zp;     
        
%% final method
            
%%from SGR to ENU
%to go back to ENU a counterclock-wise rotation is needed
pENU    = rotation_counterclockwise_mat(pSGR, TRAp.phi_ENU2local, 'rad');   
Rx_ENU  = rotation_counterclockwise(TRAp.Rxlocal, TRAp.phi_ENU2local, 'rad');   
VRx_ENU = rotation_counterclockwise(TRAp.VRxlocal, TRAp.phi_ENU2local, 'rad');  
            
%Evaluation of ENU origin in geodetic co
[ph,  lm,  ~] = xyz2llh(TRAp.OENUinECEF(1), TRAp.OENUinECEF(2), TRAp.OENUinECEF(3));
            
%Compute rotation matrix: from ecef 2 enu
Amat_ECEF2ENU = [ -sin(lm)                cos(lm)       0;
                  -sin(ph)*cos(lm)  -sin(ph)*sin(lm)  cos(ph);...
                   cos(ph)*cos(lm)   cos(ph)*sin(lm)  sin(ph)];
            
vdiff         =  [0,0,0].' -  TRAp.OENUinECEF;
OecefINenu    = Amat_ECEF2ENU*vdiff;
            
Amat_ENU2ECEF = Amat_ECEF2ENU.';
            
pdiffENU(:,:,1) = pENU(:,:,1) - OecefINenu(1);
pdiffENU(:,:,2) = pENU(:,:,2) - OecefINenu(2);
pdiffENU(:,:,3) = pENU(:,:,3) - OecefINenu(3);

pRxdiffENU      = Rx_ENU - OecefINenu;
         
%da ENU a EFEF - Evaluate to test the conversion
%pECEF   = Amat_ENU2ECEF*pdiffENU;

pECEF(:,:,1) = -sin(lm).*pdiffENU(:,:,1) - sin(ph)*cos(lm).*pdiffENU(:,:,2) + cos(ph)*cos(lm).*pdiffENU(:,:,3);
pECEF(:,:,2) =  cos(lm).*pdiffENU(:,:,1) - sin(ph)*sin(lm).*pdiffENU(:,:,2) + cos(ph)*sin(lm).*pdiffENU(:,:,3);
pECEF(:,:,3) =            0              + cos(ph).*pdiffENU(:,:,2)         + sin(ph).*pdiffENU(:,:,3);  

RxECEF  = Amat_ENU2ECEF*pRxdiffENU;
vRxECEF = Amat_ENU2ECEF*VRx_ENU;

%Conversion test
if round(RxECEF,4)  == round(TRAp.RxECEF,4),   else, disp('rf transform error'), return, end
if round(vRxECEF,4) == round(TRAp.VRxECEF,4),  else, disp('rf transform error'), return, end

%-----------------------------------------------------------------%        
%Generation of the reference system matrices to go from ECEF to OF

% TRAp.RxAttitude_angles for platfrom , TRAp.RxEuler_angles for antenna


if geoSYSp.Rx_Attitude_Variation_RPY(1) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(2) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(3) == 1
RxAttitude = [TRAp.Rx_roll,TRAp.Rx_pitch,0];
else
RxAttitude =  TRAp.Rx_nomAtt_angles;
end
[Amat_ECEF2OF, Amat_OF2BF, Amat_BFtoAF] = matrixFrameGen(RxECEF, vRxECEF, RxAttitude, TRAp.RxAttitude_angles, TRAp.RxEuler_angles);
        
%diff vector
%rs_ECEF = (pECEF - RxECEF);

%rs_OF = Amat_ECEF2OF*rs_ECEF;
%rs_BF = Amat_OF2BF*rs_OF;
%rs_AF = Amat_BFtoAF*rs_BF;            

rs_ECEF(:,:,1) = (pECEF(:,:,1) - RxECEF(1));
rs_ECEF(:,:,2) = (pECEF(:,:,2) - RxECEF(2));
rs_ECEF(:,:,3) = (pECEF(:,:,3) - RxECEF(3));
        
rs_OF(:,:,1) = Amat_ECEF2OF(1,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(1,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(1,3).*rs_ECEF(:,:,3);
rs_OF(:,:,2) = Amat_ECEF2OF(2,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(2,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(2,3).*rs_ECEF(:,:,3);
rs_OF(:,:,3) = Amat_ECEF2OF(3,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(3,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(3,3).*rs_ECEF(:,:,3);

rs_BF(:,:,1) = Amat_OF2BF(1,1).*rs_OF(:,:,1) + Amat_OF2BF(1,2).*rs_OF(:,:,2) + Amat_OF2BF(1,3).*rs_OF(:,:,3);
rs_BF(:,:,2) = Amat_OF2BF(2,1).*rs_OF(:,:,1) + Amat_OF2BF(2,2).*rs_OF(:,:,2) + Amat_OF2BF(2,3).*rs_OF(:,:,3);
rs_BF(:,:,3) = Amat_OF2BF(3,1).*rs_OF(:,:,1) + Amat_OF2BF(3,2).*rs_OF(:,:,2) + Amat_OF2BF(3,3).*rs_OF(:,:,3);

rs_AF(:,:,1) = Amat_BFtoAF(1,1).*rs_BF(:,:,1) + Amat_BFtoAF(1,2).*rs_BF(:,:,2) + Amat_BFtoAF(1,3).*rs_BF(:,:,3);
rs_AF(:,:,2) = Amat_BFtoAF(2,1).*rs_BF(:,:,1) + Amat_BFtoAF(2,2).*rs_BF(:,:,2) + Amat_BFtoAF(2,3).*rs_BF(:,:,3);
rs_AF(:,:,3) = Amat_BFtoAF(3,1).*rs_BF(:,:,1) + Amat_BFtoAF(3,2).*rs_BF(:,:,2) + Amat_BFtoAF(3,3).*rs_BF(:,:,3);
           
%% Piercing Angles in degrees
        
%Elevation angle
thetas_el_pier   = asind(rs_AF(:,:,3)./sqrt(dot(rs_AF,rs_AF,3)));
%Zenithal angle: Dot product with z axis
%thetas_zen_pier  = acosd(dot(rs_AF/norm(rs_AF), [0,0,1]));
%A cos ratio
thetas_zen_pier  = acosd(rs_AF(:,:,3)./sqrt(dot(rs_AF,rs_AF,3))); 
        
%Testing the evaluation comparing the two different definitions
%if round(thetas_zen_pier2,14) == round(thetas_zen_pier2,14), else, disp('zen piercing angles error'), return, end
        
%Azimuth angle
phis_pier       = atan2d(rs_AF(:,:,2), rs_AF(:,:,1));
                
%retrieving the index of the piercing angles
%[~, phi_ind] = (min(abs(TRAp.phiRxDown    - phis_pier)));
%[~, the_ind] = (min(abs(TRAp.thetaRxDown  - thetas_zen_pier))); 

%Gain_onsurf_lin(indy,indx)    = TRAp.GainRxDown(phi_ind, the_ind);

phiRxDown_3D      = reshape(TRAp.phiRxDown,1,1,length(TRAp.phiRxDown));
phiRxDown_array   = repmat(phiRxDown_3D,Ny_surf,Nx_surf,1);
          
thetaRxDown_3D    = reshape(TRAp.thetaRxDown,1,1,length(TRAp.thetaRxDown));
thetaRxDown_array = repmat(thetaRxDown_3D,Ny_surf,Nx_surf,1);

[~, phi_ind] = min(abs(phiRxDown_array    - phis_pier),[],3);
[~, the_ind] = min(abs(thetaRxDown_array  - thetas_zen_pier),[],3);        
indG1D       = sub2ind([size(TRAp.GainRxDown_LHCP,1), size(TRAp.GainRxDown_LHCP,2)], reshape(phi_ind, Nx_surf*Ny_surf, 1), reshape(the_ind, Nx_surf*Ny_surf, 1)); 

% %%
% 
% linear_idx = sub2ind(size(TRAp.GainRxDown_LHCP), phi_ind, the_ind);
% gainrtrt = reshape(TRAp.GainRxDown_LHCP(linear_idx), size(phi_ind));
% 
% gainrtrt = TRAp.GainRxDown_LHCP(phi_ind,the_ind);
%%

for indexTxFrequency=1:size(TRAp.GainRxDown_LHCP,3)
    Gain_LHCP_1D         = reshape(TRAp.GainRxDown_LHCP(:,:,indexTxFrequency), size(TRAp.GainRxDown_LHCP(:,:,indexTxFrequency),1).*size(TRAp.GainRxDown_LHCP(:,:,indexTxFrequency),2),1);
    Gain_LHCP_onsurf1D   = Gain_LHCP_1D(indG1D);
    Gain_LHCP_onsurf_lin = reshape(Gain_LHCP_onsurf1D, Ny_surf, Nx_surf);
    DDMparam.Gain_LHCP_onsurf_lin(:,:,indexTxFrequency) = Gain_LHCP_onsurf_lin;

    Gain_RHCP_1D         = reshape(TRAp.GainRxDown_RHCP(:,:,indexTxFrequency), size(TRAp.GainRxDown_RHCP(:,:,indexTxFrequency),1).*size(TRAp.GainRxDown_RHCP(:,:,indexTxFrequency),2),1);
    Gain_RHCP_onsurf1D   = Gain_RHCP_1D(indG1D);
    Gain_RHCP_onsurf_lin = reshape(Gain_RHCP_onsurf1D, Ny_surf, Nx_surf);
    DDMparam.Gain_RHCP_onsurf_lin(:,:,indexTxFrequency) = Gain_RHCP_onsurf_lin;

    %% Cross-pol Gain
    Gain_LR_1D         = reshape(TRAp.GainRxDown_LR(:,:,indexTxFrequency), size(TRAp.GainRxDown_LR(:,:,indexTxFrequency),1).*size(TRAp.GainRxDown_LR(:,:,indexTxFrequency),2),1);
    Gain_LR_onsurf1D   = Gain_LR_1D(indG1D);
    Gain_LR_onsurf_lin = reshape(Gain_LR_onsurf1D, Ny_surf, Nx_surf);
    DDMparam.Gain_LR_onsurf_lin(:,:,indexTxFrequency) = Gain_LR_onsurf_lin;

    Gain_RL_1D         = reshape(TRAp.GainRxDown_RL(:,:,indexTxFrequency), size(TRAp.GainRxDown_RL(:,:,indexTxFrequency),1).*size(TRAp.GainRxDown_RL(:,:,indexTxFrequency),2),1);
    Gain_RL_onsurf1D   = Gain_RL_1D(indG1D);
    Gain_RL_onsurf_lin = reshape(Gain_RL_onsurf1D, Ny_surf, Nx_surf);
    DDMparam.Gain_RL_onsurf_lin(:,:,indexTxFrequency) = Gain_RL_onsurf_lin;


    indx_origin = Xplot == 0.0;
    indy_origin = Yplot == 0.0;

    gammaParam.Gain_LHCP_onsurf_lin_atSP(indexTxFrequency) = Gain_LHCP_onsurf_lin(indy_origin, indx_origin);
    gammaParam.Gain_RHCP_onsurf_lin_atSP(indexTxFrequency) = Gain_RHCP_onsurf_lin(indy_origin, indx_origin);
    gammaParam.Gain_LR_onsurf_lin_atSP(indexTxFrequency) = Gain_LR_onsurf_lin(indy_origin, indx_origin);
    gammaParam.Gain_RL_onsurf_lin_atSP(indexTxFrequency) = Gain_RL_onsurf_lin(indy_origin, indx_origin);

%% -------- DOPPLER SHIFT COMPUTATION ----------------------------
      
%           if TRAp.Rxlocal(3) < 15e3
%           vRxECEF_EARTH_ROT = TRAp.VRxECEF + (cross([0,0,7.2921158553*1.e-5],TRAp.RxECEF)).';
%           vRxENU_EARTH_ROT = ecef2enu_dc(vRxECEF_EARTH_ROT,TRAp.OENUinECEF, 'v');   
%           VRx_local  = rotation_clockwise(vRxENU_EARTH_ROT,   TRAp.phi_ENU2local, 'rad');
% %           vHydroG_ECEF  = vHydroG_ECEF + (cross([0,0,7.2921158553*1.e-5],pHydroG_ECEF)).';
%          end
%       DDS = transpose(RTx)*(VRx_local - TRAp.VTxlocal); % Speed diff * versor * WaveL [Hertz]
%     TRAp.VRxlocal = [0 0 0]';
    DDS = transpose(RTx)*(TRAp.VRxlocal - TRAp.VTxlocal); % Speed diff * versor * WaveL [Hertz]
    
%     DDS = transpose(RTx)*(TRAp.VRxlocal - TRAp.VTxlocal); % Speed diff * versor * WaveL [Hertz]

    DDS = DDS/RTx_mod/TRAp.waveL(indexTxFrequency);   

    DDR = refl_doppler2_GNSSR_mat(pSGR, TRAp, TRAp.waveL(indexTxFrequency));
  
    %doppler offset (to keep the SP at doppler shift = 0)
    [DDR_SP,VRx_local_earth_rot] = refl_doppler2_GNSSR(TRAp.phi_ENU2local,0.0, 0.0, 0.0, TRAp.Txlocal, TRAp.VTxlocal, TRAp.Rxlocal, TRAp.VRxlocal, TRAp.OENUinECEF, TRAp.waveL(indexTxFrequency),TRAp);

    DDR_offset = DDR_SP - DDS;
    % [Hertz] Doppler shift taking into account the offset 
    DDMparam.fdoppler(:,:,indexTxFrequency) = DDR - DDS - DDR_offset;


%% --------------- DELAY COMPUTATION-------------------------------

    %Distance Rx to point on ground
    RxPdiff(:,:,1) = TRAp.Rxlocal(1) - pSGR(:,:,1); 
    RxPdiff(:,:,2) = TRAp.Rxlocal(2) - pSGR(:,:,2); 
    RxPdiff(:,:,3) = TRAp.Rxlocal(3) - pSGR(:,:,3); 

    %norm Rx-P
    RxPmod  = sqrt(RxPdiff(:,:,1).^2 + RxPdiff(:,:,2).^2 + RxPdiff(:,:,3).^2);          
        
    %Distance Tx to point on ground
    TxPdiff(:,:,1) = TRAp.Txlocal(1) - pSGR(:,:,1); 
    TxPdiff(:,:,2) = TRAp.Txlocal(2) - pSGR(:,:,2); 
    TxPdiff(:,:,3) = TRAp.Txlocal(3) - pSGR(:,:,3); 

    %norm Tx-P    
    TxPmod  = sqrt(TxPdiff(:,:,1).^2 + TxPdiff(:,:,2).^2 + TxPdiff(:,:,3).^2);

    %same extra path due to ionosphere for all points on the ground
    DDMparam.fdelay(:,:,indexTxFrequency)  = (RxPmod + TxPmod + ionosphericDelay(indexTxFrequency) - refdelay(indexTxFrequency))./speedlight;              % [microseconds: usec] - linee isodelay
end      
  
%% -------scattering angles at each pixel point (xp,yp, zp) in rotated ENU------
%  ----- Note: theta_s is an incidence angle (zenithal), and as such
%                will be used later on in the code
       
sigmaParam.TH_s = asind(sqrt((pSGR(:,:,1) - TRAp.Rxlocal(1)).^2 + (pSGR(:,:,2) - TRAp.Rxlocal(2)).^2)./RxPmod);
sigmaParam.PH_s = atan2d(TRAp.Rxlocal(2)  - pSGR(:,:,2), TRAp.Rxlocal(1) - pSGR(:,:,1));  

% %alternative - it is a matrix - not practical to test
% zax(:,:,1) = zeros(Ny_surf,Nx_surf);
% zax(:,:,2) = zeros(Ny_surf,Nx_surf);
% zax(:,:,3) = ones(Ny_surf,Nx_surf);
% 
% TH_S2  = acosd(dot(RxPdiff./RxPmod, zax, 3));
% TH0_S2 = acosd(dot(TxPdiff./TxPmod, zax, 3));        
        

%% ********Incident and scattered polarization unit vectors in rotated ENU
        
[pt, prR, prL] = polarization_unit_vectors_mat(TRAp, sigmaParam.TH_s, 'deg');      
        

%% Implemented for Rongowai and GREAT - AMIR
% Instead of calculating the polarization unit vector, for the input of SAVERS it
% has been considered an ideal antenna. Based on the cross-polarization pattern, it
% will be considered in the calculation of power.
pt(1) = 1;
pt(2) = complex(0,-1);
prR(:,:,1) = ones(Ny_surf, Nx_surf); 
prR(:,:,2) = complex(0,-1).*ones(Ny_surf, Nx_surf);
prL(:,:,1) = ones(Ny_surf, Nx_surf);
prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);
%%

sigmaParam.ptreal    = real(pt);
sigmaParam.ptim      = imag(pt);
        
%%%%Attenzione qui, la terza dimensione era inutilizzata secondo me
sigmaParam.prR_real(:,:,1) = real(prR(:,:,1));
sigmaParam.prR_real(:,:,2) = zeros(Ny_surf,Nx_surf);

sigmaParam.prR_im(:,:,1)   = zeros(Ny_surf,Nx_surf);
sigmaParam.prR_im(:,:,2)   = imag(prR(:,:,2));

sigmaParam.prL_real(:,:,1) = real(prL(:,:,1));
sigmaParam.prL_real(:,:,2) = zeros(Ny_surf,Nx_surf);

sigmaParam.prL_im(:,:,1)   = zeros(Ny_surf,Nx_surf);
sigmaParam.prL_im(:,:,2)   = imag(prL(:,:,2));
       

%%   incidence and scattering angles in the facet local reference (FLR)
        
       
%incidence and scattering angles in the facet local reference (FLR) at each
%pixel point (xp,yp,zp), that is the ENU rotated frame which has been
%also tilted to take into account of the slope and aspect of the facet
        
alpha               = 90.0 - aspect;
alpha(isnan(alpha)) = 0.0;
beta                = slope;
beta(isnan(beta))   = 0.0;
        
%[incAngle_flr(indy,indx),incAzimuthAngle_flr(indy,indx),TH_s_flr(indy,indx),PH_s_flr(indy,indx)] = ...
%                    transformation_abs_flr(beta, alpha, TRAp.incAngle_enuRot, TH_s(indy,indx), PH_s(indy,indx));

[sigmaParam.incAngle_flr, sigmaParam.incAzimuthAngle_flr, sigmaParam.TH_s_flr, sigmaParam.PH_s_flr] = transformation_abs_flr_mat(...
                                                                                                              beta, alpha, TRAp.incAngle_enuRot,...
                                                                                                              sigmaParam.TH_s, sigmaParam.PH_s);        
        
%Unit vector of the normal at (xp, yp, zp)
%zVect_flr(indy,indx,:)  = [sind(beta)*sind(alpha); -cosd(alpha)*sind(beta); cosd(beta)];
sigmaParam.zVect_flr(:,:,1)  =  sind(beta).*sind(alpha);
sigmaParam.zVect_flr(:,:,2)  = -cosd(alpha).*sind(beta); cosd(beta);
sigmaParam.zVect_flr(:,:,3)  =  cosd(beta);

%deltaPhi_flr(indy,indx) =  wrapTo180(PH_s_flr(indy,indx)-incAzimuthAngle_flr(indy,indx));
sigmaParam.deltaPhi_flr =  wrapTo180(sigmaParam.PH_s_flr - sigmaParam.incAzimuthAngle_flr);
        
%Computing path lenght weighting factor at each point (xp,yp)
DDMparam.PATH_factor = 1./(RxPmod.*TxPmod);    % Bistatic Path Factor
DDMparam.TxPmod      = TxPmod;
DDMparam.RxPmod      = RxPmod;
sigmaParam.r2mat     = RxPmod;              % distance vs.Rx [m]     

% PATH_factor_SP = (RxPmod.*TxPmod); %AMIR test
% gammaParam.PATH_factor_SP =  PATH_factor_SP(ceil(Nx_surf/2),ceil(Ny_surf/2));     
%close(waitbarSimGrid) %======close waitbar

%------------------------------------------------------------------------

%% plot theta, phi,gain, isolinee

%% No error - No plot for error free part of simulation.
if ~isfield(TRAp, 'InternalFlag_NoError') || TRAp.InternalFlag_NoError ~= 1

if tag.plotAngles_Isolines == 1

    %figure isolines
    figIsoLinesPattern = figure('Name', 'Isolines L1', 'NumberTitle','off');
    %set(gcf,'Position',[0 0 750 750])

    subplot(2,2,1) % upper left
    imagesc(Xplot,Yplot,sigmaParam.incAngle_flr)
    axis xy
    cbarDEM=colorbar;
    cbarDEM.Label.String = 'Incidence angle  (deg)';
    title('Local incidence angle')
    xlabel('X axis (m)')
    ylabel('Y axis (m)')   

    subplot(2,2,2) % upper right
    imagesc(Xplot,Yplot,sigmaParam.TH_s_flr)
    axis xy
    cbarDEM=colorbar;
    cbarDEM.Label.String = 'angle  (deg)';
    title('Local scattering angle')
    xlabel('X axis (m)')
    ylabel('Y axis (m)')  

    subplot(2,2,3) % lower left
    [~,h] = contour(Xplot,Yplot, DDMparam.Gain_LHCP_onsurf_lin(:,:,1)/max(max(DDMparam.Gain_LHCP_onsurf_lin(:,:,1))));
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
    axis('equal')
    title('Rx antenna pattern for L1 LHCP')
    xlabel('X axis [m]')
    ylabel('Y axis [m]')        

    subplot(2,2,4) %lower right
    hold
    [~,h] = contour(Xplot,Yplot,DDMparam.fdoppler(:,:,1));% <<<<<<<<<<<<<<<<< [Hertz]
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
    [~,h] = contour(Xplot,Yplot,DDMparam.fdelay(:,:,1));  % <<<<<<<<<<<<<<<<< [microsec]
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
    hold
    axis('equal')
    title('        IsoDop.(Hz) IsoRan.(usec)')
    xlabel('X axis [m]'), ylabel('Y axis [m]')

    saveas(figIsoLinesPattern, strjoin([nameOutputFile '_isoLinesPatterns_L1.fig'],'')); 

    %Amir
    figIsolines_antFootprint = figure( ...
    'Name', 'Isolines Vs. Antenna Footprint L1', ...
    'NumberTitle', 'off');
    hold on;

    [C1, h1] = contour(Xplot, Yplot, DDMparam.fdoppler(:,:,1), 'r');
%     set(h1, 'ShowText', 'on');  % add numeric labels to contours
    Zdelay = DDMparam.fdelay(:,:,1);
    levels_delay = 1 : ceil(max(Zdelay(:)));
    [C2, h2] = contour(Xplot, Yplot, Zdelay, levels_delay, 'b');
%      set(h2, 'ShowText', 'on');  % add numeric labels to contours
    [~, h2_bold] = contour(Xplot, Yplot, Zdelay, [1 1], 'b', 'LineWidth', 2.5);
    [C_delay, ~] = contour(Xplot, Yplot, Zdelay, [1 1], 'b', 'LineWidth', 2.5);

    % Calculate first chip zone semi major axis

    contour_points_delay = [];
    k = 1;
    while k < size(C_delay, 2)
        num_points = C_delay(2, k);
        x = C_delay(1, k+1:k+num_points);
        y = C_delay(2, k+1:k+num_points);
        contour_points_delay = [contour_points_delay; x', y'];
        k = k + num_points + 1;
    end

    D_delay = pdist(contour_points_delay);  
    diameter_delay_max = max(D_delay);


%     levels_to_show = 0.5:0.1:1.0;
%     [C3, h3] = contour(Xplot, Yplot, gain_norm, levels_to_show, 'k');

     % Gain contours
    gain_norm = DDMparam.Gain_LHCP_onsurf_lin(:,:,1) / max(max(DDMparam.Gain_LHCP_onsurf_lin(:,:,1)));
    
    max_gain_in_plot = 0.5;
    [C_bold, h_bold] = contour(Xplot, Yplot, gain_norm,[max_gain_in_plot max_gain_in_plot], 'k', 'LineWidth', 3);
    
    if isempty(C_bold) == 1
        for max_gain_in_plot = 0.6:0.1:0.9
            [C_bold, h_bold] = contour(Xplot, Yplot, gain_norm, [max_gain_in_plot max_gain_in_plot], 'k', 'LineWidth', 3);
            
            if isempty(C_bold) == 0
                 break; 
            end
        end
    end

    min_visible_gain = ceil(min(gain_norm(:)) * 10) / 10;
    [C_dashed, h_dashed] = contour(Xplot, Yplot, gain_norm, [min_visible_gain min_visible_gain], ...
                                   'k', 'LineWidth', 2);
    set(h_dashed, 'LineStyle', '--');
    clabel(C_dashed, h_dashed, 'Color', 'k', 'FontSize', 15);

    % Calculate diameter gain norm 0.5 contour
    contour_points = [];
    k = 1;
    while k < size(C_bold, 2)
        num_points = C_bold(2, k);
        x = C_bold(1, k+1:k+num_points);
        y = C_bold(2, k+1:k+num_points);
        contour_points = [contour_points; x', y'];
        k = k + num_points + 1;
    end

    D = pdist(contour_points);
    diameter_max = max(D);

    p1 = patch(NaN, NaN, 'r', 'DisplayName', 'IsoDoppler');
    p2 = patch(NaN, NaN, 'b', 'DisplayName', 'IsoDelay');
    p3 = plot(NaN, NaN, 'k-', 'LineWidth', 3, ...
          'DisplayName', sprintf('Rx antenna footprint (%.1f power)', max_gain_in_plot));
    p4 = plot(NaN, NaN, 'k--', 'LineWidth', 3, 'DisplayName', 'Rx beam low gain edge');

    legend([p1 p2 p3 p4], 'Location', 'southwest');

    %% AMIR TEST HEADING

Rx_xy = TRAp.Rxlocal(1:2);          % take only x and y components

%draw the receiver as a filled black dot
p5 = plot(Rx_xy(1), Rx_xy(2), 'ko', ...
          'MarkerSize', 6, 'MarkerFaceColor', 'k', ...
          'DisplayName', 'Rx position');

% Rx_heading_Local_rad_earth_rot = atan2(VRx_local_earth_rot(2), VRx_local_earth_rot(1));
% u =  cos(Rx_heading_Local_rad_earth_rot);   % X local component
% v =  sin(Rx_heading_Local_rad_earth_rot);   % Y local component

% heading to E,N
heading_rad = deg2rad(TRAp.Rx_heading_Local_deg);
u =  cos(heading_rad);   % X local component
v =  sin(heading_rad);   % Y local component



if isempty(diameter_max) == 1
    diameter_max = Xplot(end)/2;
end
arrow_len = 0.3 * diameter_max;    % % of the HP-footprint diameter

%raw the arrow with quiver; ‘0’ turns off autoscaling
p6 = quiver(Rx_xy(1), Rx_xy(2), ...
            arrow_len*u, arrow_len*v, 0, ...
            'k', 'LineWidth', 2, 'MaxHeadSize', 0.7, ...
            'DisplayName', 'Rx heading');

% %heading to E,N
% heading_rad_gps = deg2rad(TRAp.Tx_heading_Local_deg);
% u_tx =  cos(heading_rad_gps);   % X local component
% v_tx =  sin(heading_rad_gps);   % Y local component
% 
% p7 = quiver(Rx_xy(1), Rx_xy(2), ...
%             arrow_len*u_tx, arrow_len*v_tx, 0, ...
%             'Color', 'm', 'LineWidth', 2, 'LineStyle', '-', ...
%             'MaxHeadSize', 0.7, 'DisplayName', 'GPS heading');


%% ──────────────────────────  update the legend  ─────────────────────────────
legend([p1 p2 p3 p4 p5 p6], 'Location', 'southwest');

    %% AMIR TEST HEADING

    xlabel('X');
    ylabel('Y');
    title('IsoDoppler, IsoDelay, and Antenna Footprint L1');
    subtitle(sprintf( ...
    'Diameter Half Power ANT Footprint: %.2f m\nMajor axis First Chip Zone: %.2f m', ...
    diameter_max, diameter_delay_max));
    axis equal;
    grid on;

    saveas(figIsolines_antFootprint, strjoin([nameOutputFile '_isoLines_Vs_AntFoortprint_L1.fig'],'')); 

    %Amir
  %  if geoSYSp.Mode ~= "Single Point" , close (figIsoLinesPattern), end ; 

    if length(TRAp.waveL)==2
        %figure isolines
        figIsoLinesPattern_L5 = figure('Name', 'Isolines L5', 'NumberTitle','off');        
        subplot(2,2,3) % lower left
        [~,h] = contour(Xplot,Yplot, DDMparam.Gain_LHCP_onsurf_lin(:,:,2)/max(max(DDMparam.Gain_LHCP_onsurf_lin(:,:,2))));
        set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
        axis('equal')
        title('Rx antenna pattern for L5 LHCP')
        xlabel('X axis [m]')
        ylabel('Y axis [m]')        

        subplot(2,2,4) %lower right
        hold
        [~,h] = contour(Xplot,Yplot,DDMparam.fdoppler(:,:,2));% <<<<<<<<<<<<<<<<< [Hertz]
        set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
        [~,h] = contour(Xplot,Yplot,DDMparam.fdelay(:,:,2));  % <<<<<<<<<<<<<<<<< [microsec]
        set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
        hold
        axis('equal')
        title('        IsoDop.(Hz) IsoRan.(usec)')
        xlabel('X axis [m]'), ylabel('Y axis [m]')

        saveas(figIsoLinesPattern_L5, strjoin([nameOutputFile '_isoLinesPatterns_L5.fig'],''));
        close (figIsoLinesPattern_L5)

            %Amir
    figIsolines_antFootprint_L5= figure('Name', 'Isolines Vs. Antenna Footprint L5', 'NumberTitle','off');
    hold on;

    [C1, ~] = contour(Xplot, Yplot, DDMparam.fdoppler(:,:,2), 'r');

    Zdelay5 = DDMparam.fdelay(:,:,2);
    levels_delay = 1 : ceil(max(Zdelay5(:)));
    [C2, h2] = contour(Xplot, Yplot, Zdelay5, levels_delay, 'b');
    set(h2, 'ShowText', 'on');  % add numeric labels to contours
    [~, h2_bold] = contour(Xplot, Yplot, Zdelay5, [1 1], 'b', 'LineWidth', 2.5);
     gain_norm5 = DDMparam.Gain_LHCP_onsurf_lin(:,:,2) / max(max(DDMparam.Gain_LHCP_onsurf_lin(:,:,2)));

    [~, h_bold]   = contour(Xplot, Yplot, gain_norm5, [0.5 0.5], 'k',  'LineWidth', 3);     % solid thick line
    min_visible_gain5 = ceil(min(gain_norm5(:)) * 10) / 10;
    [C_dashed, h_dashed] = contour(Xplot, Yplot, gain_norm5, [min_visible_gain5 min_visible_gain5], ...
                        'k', 'LineWidth', 2);
    set(h_dashed, 'LineStyle', '--');

    clabel(C_dashed, h_dashed, 'Color', 'k', 'FontSize', 15);

    % Legend entries (use plot for lines, patch for areas)
    p1 = patch(NaN, NaN, 'r', 'DisplayName', 'IsoDoppler');
    p2 = patch(NaN, NaN, 'b', 'DisplayName', 'IsoDelay');
    p3 = plot(NaN, NaN, 'k-', 'LineWidth', 3, 'DisplayName', 'Rx antenna footprint (Half Power)');
    p4 = plot(NaN, NaN, 'k--', 'LineWidth', 3, 'DisplayName', 'Rx beam low gain edge');

    legend([p1 p2 p3 p4], 'Location', 'southwest');
    xlabel('X');
    ylabel('Y');
    title('IsoDoppler, IsoDelay, and Antenna Footprint L5');
    axis equal;
    grid on;

    saveas(figIsolines_antFootprint_L5, strjoin([nameOutputFile '_isoLines_Vs_AntFoortprint_L5.fig'],'')); 
%     close (figIsolines_antFootprint_L5)
    %Amir
    end
    
end
end
   
end