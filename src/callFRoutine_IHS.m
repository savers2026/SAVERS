%% ************************************************************************
% MODULE NAME:      callFortranRoutine_GNSSR_inhomogSurface.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Dente, D. Comite, N. Pierdicca, L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      2018 by Leila Guerriero and Roberto Giusto
% Main updates:          2020 by L. Dente
%                        2022 by D. Comite
%                        Mar 2023 by Laura Dente (L1&L5, E1&E5 in parallel)
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module calls the Fortran routines that computes the sigma0 according
% to the AIEM for the bare soils, the forest model by Tor Vergata and the
% low vegetation model by Tor Vergata. It computes the polarization loss
% factor fot he case of vegetated surfaces (note that for the case of bare
% soils this is done inside the Fortran routine).
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [sigmaRR_cohe,sigmaLR_cohe,sigmaRR_inco,sigmaLR_inco] = callFRoutine_IHS(Nx_surf,Ny_surf,...
                                                                 TxBeamwidth,TxFrequency,RxBeamwidth,incAngle_enuRot,sigmaParam,...
                                                                 terrainParameters_enuRot,vegetationParameters_enuRot,TRAp,gammaParam,surfParam,geoSYSp)

beta0rad =  deg2rad(TxBeamwidth);
beta2rad =  deg2rad(RxBeamwidth);

terrainFlag        = terrainParameters_enuRot(1);
sds_ini            = terrainParameters_enuRot(2);
l_ini              = terrainParameters_enuRot(3);
corr_funct         = terrainParameters_enuRot(4);
mv_ini             = terrainParameters_enuRot(5);
bulkDensity        = terrainParameters_enuRot(6);
surfaceTemperature = terrainParameters_enuRot(7);
                                        
indexPlant   = vegetationParameters_enuRot(1);
kindPlant    = vegetationParameters_enuRot(2);
LAI          = vegetationParameters_enuRot(3);
amoipl       = vegetationParameters_enuRot(4);
wleaf        = vegetationParameters_enuRot(5);
bleaf        = vegetationParameters_enuRot(6);
Bio          = vegetationParameters_enuRot(7);
NTree        = vegetationParameters_enuRot(8);
dbh_mean     = vegetationParameters_enuRot(9);
DbhStD       = vegetationParameters_enuRot(10);
kindDBH_dist = vegetationParameters_enuRot(11);
allomEq      = vegetationParameters_enuRot(12);


% sds_ini = 0;
% AMIR TEST Betac
if RxBeamwidth == 0

[~ , ~ , R0] = xyz2llh(TRAp.TxECEF(1),TRAp.TxECEF(2),TRAp.TxECEF(3));
[~ , ~ , Rs] = xyz2llh(TRAp.RxECEF(1),TRAp.RxECEF(2),TRAp.RxECEF(3));


% beta_opt = findBeta0Beta2_theta(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam,geoSYSp);

% beta_opt = findBeta0Beta2_New(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam);
% beta_opt = findBeta0Beta2_theta_i(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam);
beta_opt = findBeta0Beta2_performance(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam,geoSYSp);
% beta_opt = findBeta0Beta2_Rs(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam,geoSYSp);

% beta_opt = findBeta0Beta2_New_Coarse(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam);
% beta_c = findBetaC_New(R0, Rs, TRAp.waveL, TRAp.TxSignalPower, sds_ini, TRAp.GainRxDown_LHCP  , Ny_surf, Nx_surf, TRAp.GainTx_linear, sigmaParam, TRAp, mv_ini,gammaParam,surfParam);
beta2rad =  deg2rad(beta_opt(2));
beta0rad =  deg2rad(beta_opt(1));

end
% AMIR TEST Betac

%% --- SIGMA0 COMPUTATION for not vegetated surfaces ---
%run the bare soil routine for the whole integration area
%it is temporary assumed that the complete integration area has a cover
%class equal to the coverageID given as input to the routine (as well as
%all the soil and vegetation parameters defined for that class)
%This step is necessary because of the resampling done on TH_S and PH_s
%In the routine SAVERS_spaceborneGNSSR_inhomogSurface only the facets with the bare soil flag on will be considered
%and a mosaic of the sigma0 will be carried out

if indexPlant == 0 || LAI==0 && Bio==0


    [sigmaRR_cohe,sigmaRR_inco,sigmaLR_cohe,sigmaLR_inco] = AIEM_GNSSRFLR_DELPHIS_mironov_covariance_IF(Nx_surf,Ny_surf,...
                                                            incAngle_enuRot,sigmaParam.TH_s,sigmaParam.PH_s,sigmaParam.zVect_flr,...
                                                            sigmaParam.incAngle_flr,sigmaParam.TH_s_flrIF,sigmaParam.deltaPhi_flr,...
                                                            sigmaParam.ptreal,sigmaParam.ptim,sigmaParam.prR_real,sigmaParam.prR_im,...
                                                            sigmaParam.prL_real,sigmaParam.prL_im,...
                                                            mv_ini,bulkDensity,sds_ini,l_ini,corr_funct,terrainFlag,surfaceTemperature,...
                                                            beta0rad,beta2rad,sigmaParam.r2mat,sigmaParam.R0,TxFrequency);

else

%% --- SIGMA ZERO COMPUTATION FOR VEGETATED SURFACES ---
%run the vegetation routine for the whole integration area
%it is temporary assumed that the complete integration area has a cover
%class equal to the coverageID given as input to the routine (as well as
%all the soil and vegetation parameters defined for that class)
%This step is necessary because of the resampling done on TH_S and PH_s
%In the routine combineCoverageContributions only the facets with the vegetation flag on will be considered

    if kindPlant == 1

        %define descrete angles to compute the volume scattering
        % (i) row index from 5° to 85°  Theta_angle of scattering
        Theta_nx_forest = 18; %FOREST
        % (j) column index from 0° to 180°  Phi_angle of scattering
        Phi_ny_forest = 33; %FOREST
        TH_srfor      = zeros(Phi_ny_forest,Theta_nx_forest);    % reduced matrix scattering in Elevation
        PH_srfor      = zeros(Phi_ny_forest,Theta_nx_forest);    % reduced matrix scattering in Azimuth
              
        if allomEq == 0 %broadleaved forest
            [sigmaRR_inco_noDepol,sigmaRR_cohe_noDepol,ForsRR_inco,...
             sigmaLR_inco_noDepol,sigmaLR_cohe_noDepol,ForsLR_inco] = Forest_new_spaceborneGNSSR_vegloc_mironov(Nx_surf,Ny_surf,...
                                                                      incAngle_enuRot,sigmaParam.TH_s,sigmaParam.PH_s,sigmaParam.zVect_flr,...
                                                                      sigmaParam.incAngle_flr, sigmaParam.TH_s_flr,sigmaParam.deltaPhi_flr,...
                                                                      mv_ini,bulkDensity,sds_ini,l_ini,terrainFlag,surfaceTemperature,Bio,kindDBH_dist,LAI,dbh_mean,DbhStD,...
                                                                      NTree,beta0rad,beta2rad,sigmaParam.r2mat,sigmaParam.R0,TxFrequency);
            clear Forest_new_spaceborneGNSSR_TDS_vegloc_mironov
        elseif allomEq==1 %not used in HSAVERS
            [sigmaRR_inco_noDepol,sigmaRR_cohe_noDepol,ForsRR_inco,...
             sigmaLR_inco_noDepol,sigmaLR_cohe_noDepol,ForsLR_inco] = Forest_new_spaceborneGNSSR_vegloc_amaz_mironov(Nx_surf,Ny_surf,...
                                                                      incAngle_enuRot,sigmaParam.TH_s,sigmaParam.PH_s,sigmaParam.zVect_flr,...
                                                                      sigmaParam.incAngle_flr, sigmaParam.TH_s_flr,sigmaParam.deltaPhi_flr,...
                                                                      mv_ini,bulkDensity,sds_ini,l_ini,terrainFlag,surfaceTemperature,Bio,kindDBH_dist,LAI,dbh_mean,DbhStD,...
                                                                      NTree,beta0rad,beta2rad,sigmaParam.r2mat,sigmaParam.R0,TxFrequency);
            clear Forest_new_spaceborneGNSSR_TDS_vegloc_amaz_mironov
        else % needleleaved forest
            [sigmaRR_inco_noDepol,sigmaRR_cohe_noDepol,ForsRR_inco,...
             sigmaLR_inco_noDepol,sigmaLR_cohe_noDepol,ForsLR_inco] = Forest_new_spaceborneGNSSR_vegloc_softwood_mironov(Nx_surf,Ny_surf,...
                                                                      incAngle_enuRot,sigmaParam.TH_s,sigmaParam.PH_s,sigmaParam.zVect_flr,...
                                                                      sigmaParam.incAngle_flr, sigmaParam.TH_s_flr,sigmaParam.deltaPhi_flr,...
                                                                      mv_ini,bulkDensity,sds_ini,l_ini,terrainFlag,surfaceTemperature,Bio,kindDBH_dist,LAI,dbh_mean,DbhStD,...
                                                                      NTree,beta0rad,beta2rad,sigmaParam.r2mat,sigmaParam.R0,TxFrequency);
            clear Forest_new_spaceborneGNSSR_TDS_vegloc_softwood_mironov         
        end

        %interpolate the volume scattering over all angles of scattering
        [TH_srfor, PH_srfor] = meshgrid(linspace(2.5,87.5,18), linspace(0,180,33));

%         sigmaRR_veg = interp2(TH_srfor, PH_srfor, ForsRR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'spline');
%         sigmaLR_veg = interp2(TH_srfor, PH_srfor, ForsLR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'spline');

%         sigmaRR_veg = interp2(TH_srfor, PH_srfor, ForsRR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'makima');
%         sigmaLR_veg = interp2(TH_srfor, PH_srfor, ForsLR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'makima');
%         sigmaRR_veg(sigmaRR_veg<0)=0;
%         sigmaLR_veg(sigmaLR_veg<0)=0;
        sigmaRR_veg = interp2(TH_srfor, PH_srfor, ForsRR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'linear');
        sigmaLR_veg = interp2(TH_srfor, PH_srfor, ForsLR_inco', sigmaParam.TH_s, abs(sigmaParam.PH_s),'linear');

        %sum the soil incohent component to the vegetation volume component
        %(sigmaLR_cohe_noDepol is the coehent component of soil and vegetation)
        sigmaRR_inco_noDepol = sigmaRR_inco_noDepol + sigmaRR_veg;
        sigmaLR_inco_noDepol = sigmaLR_inco_noDepol + sigmaLR_veg;

        pRideal_hat(1) = 1;     % ideal Right pol.vector
        pRideal_hat(2) = complex(0,-1);
        pRideal_hat    = pRideal_hat./norm(pRideal_hat);
        pLideal_hat(1) = 1;     % ideal Left pol.vector
        pLideal_hat(2) = complex(0,1);
        pLideal_hat    = pLideal_hat./norm(pLideal_hat);

        %polarization loss factor (PLF) of co and inco components
            pRideal_hat2(:,:,1) = ones(Ny_surf, Nx_surf).*sqrt(2)./2; 
            pRideal_hat2(:,:,2) = -1j*ones(Ny_surf, Nx_surf).*sqrt(2)./2;

            pLideal_hat2(:,:,1) = ones(Ny_surf, Nx_surf).*sqrt(2)./2; 
            pLideal_hat2(:,:,2) = 1j*ones(Ny_surf, Nx_surf).*sqrt(2)./2; 

            prR_hat2(:,:,1)     = complex(sigmaParam.prR_real(:,:,1), sigmaParam.prR_im(:,:,1));
            prR_hat2(:,:,2)     = complex(sigmaParam.prR_real(:,:,2), sigmaParam.prR_im(:,:,2));
            normprR_hat2        = sqrt(abs(prR_hat2(:,:,1)).^2 + abs(prR_hat2(:,:,2)).^2);
            prR_hat2            = prR_hat2./normprR_hat2;

            prL_hat2(:,:,1)     = complex(sigmaParam.prL_real(:,:,1), sigmaParam.prL_im(:,:,1));
            prL_hat2(:,:,2)     = complex(sigmaParam.prL_real(:,:,2), sigmaParam.prL_im(:,:,2));
            normprL_hat2        = sqrt(abs(prL_hat2(:,:,1)).^2 + abs(prL_hat2(:,:,2)).^2);
            prL_hat2            = prL_hat2./normprL_hat2;

            PLF_RR = dot(prR_hat2,pRideal_hat2,3).*conj(dot(prR_hat2, pRideal_hat2,3));
            PLF_LR = dot(prL_hat2,pRideal_hat2,3).*conj(dot(prL_hat2, pRideal_hat2,3));
            PLF_RL = dot(prR_hat2,pLideal_hat2,3).*conj(dot(prR_hat2, pLideal_hat2,3));
            PLF_LL = dot(prL_hat2,pLideal_hat2,3).*conj(dot(prL_hat2, pLideal_hat2,3));      

            sigmaRR_inco = sigmaRR_inco_noDepol.*PLF_RR + sigmaLR_inco_noDepol.*PLF_RL; 
            sigmaLR_inco = sigmaLR_inco_noDepol.*PLF_LL + sigmaRR_inco_noDepol.*PLF_LR; 
            sigmaRR_cohe = sigmaRR_cohe_noDepol.*PLF_RR + sigmaLR_cohe_noDepol.*PLF_RL; 
            sigmaLR_cohe = sigmaLR_cohe_noDepol.*PLF_LL + sigmaRR_cohe_noDepol.*PLF_LR;


    else
%% --------   LOW VEGETATION     ---------------------------

        %density of leaves
        %LAI/leaf area
        an0leaf=LAI/(pi*(wleaf/2)*(bleaf/2));
        
        [sigmaRR_inco,sigmaRR_cohe,sigmaLR_inco,sigmaLR_cohe] = shrubs_covariance_dem_all_mironov(Nx_surf,Ny_surf,...
                                                                incAngle_enuRot,sigmaParam.TH_s,sigmaParam.PH_s,sigmaParam.zVect_flr,...
                                                                sigmaParam.incAngle_flr,sigmaParam.TH_s_flr,sigmaParam.deltaPhi_flr,...
                                                                sigmaParam.ptreal,sigmaParam.ptim,sigmaParam.prR_real,...
                                                                sigmaParam.prR_im,sigmaParam.prL_real,sigmaParam.prL_im,...
                                                                mv_ini,bulkDensity,sds_ini,l_ini,terrainFlag,surfaceTemperature, ...
                                                                an0leaf,amoipl,wleaf,bleaf,beta0rad,beta2rad,sigmaParam.r2mat,sigmaParam.R0,TxFrequency);
        clear shrubs_covariance_dem_all_mironov
        
%          close(hwt2) %===========================================================
    end

end
%% 
% if the incidence and/or scattering angle in the facet reference system
% are larger than 70 deg, then the sigma0 is set to NaN.

sigmaRR_cohe(sigmaParam.incAngle_flr>70.0 | sigmaParam.TH_s_flr>70.0) = NaN;
sigmaLR_cohe(sigmaParam.incAngle_flr>70.0 | sigmaParam.TH_s_flr>70.0) = NaN;
sigmaRR_inco(sigmaParam.incAngle_flr>70.0 | sigmaParam.TH_s_flr>70.0) = NaN;
sigmaLR_inco(sigmaParam.incAngle_flr>70.0 | sigmaParam.TH_s_flr>70.0) = NaN;


end


