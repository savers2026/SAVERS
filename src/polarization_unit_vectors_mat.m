%% ************************************************************************
% MODULE NAME:      polarization_unit_vectors_mat.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      D. Comite and L. Guerriero
% @copyright:   Sapienza University of Rome
% Original version:      2021-2022 by L. Guerriero
% Main updates:          2021 - 2022 by D. Comite
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% 
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

%% Name: transformation_abs_flr
%% Authors: Leila Guerriero, Tor Vergata University, Davide Comite, Sapienza University
%% New version: Mar 2022

function [pt, prR, prL] = polarization_unit_vectors_mat(TRAp, thetaRx_pier, tag)

          flagTxMismatch = TRAp.flagTxMismatch;
          flagRxmismatch = TRAp.flagRxDOWNmismatch;
          thetaTx_pier   = TRAp.incAngle_enuRot; 

          Ny_surf = size(thetaRx_pier,2);
          Nx_surf = size(thetaRx_pier,1);

switch tag
    case 'rad'
        
        if flagTxMismatch == 0
            pt(1) = 1;                                           %Transmitted ideal Right pol.vector
            pt(2) = complex(0,-1);
        else
            pt(1) = cos(thetaTx_pier);                           %Right pol.vector with mismatch
            pt(2) = complex(0,-1);
        end
        
        if flagRxmismatch == 0
            prR(:,:,1) = ones(Ny_surf, Nx_surf);                     %Received ideal Right pol.vector
            prR(:,:,2) = complex(0,-1).*ones(Ny_surf, Nx_surf);
            prL(:,:,1) = ones(Ny_surf, Nx_surf);                     %Received ideal Left pol.vector
            prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);
        else
%             prR(:,:,1) = cos(thetaRx_pier);                          %Right pol.vector with mismatch
%             prR(:,:,2) = complex(0,-1).*ones(Ny_surf, Nx_surf);
%             prL(:,:,1) = cos(thetaRx_pier);                          %Left pol.vector with mismatch
%             prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);            
            thetaRx_pier=rad2deg(thetaRx_pier);
            stdtheta= TRAp.MM_RHCP(1); bias=TRAp.MM_RHCP(2) ; %%%theta in degrees%%%
            prR(:,:,1)=funcpol(thetaRx_pier,stdtheta^2,bias);  %exp MM Right pol antenna
            prR(:,:,2)=complex(0,-1).*ones(Ny_surf, Nx_surf); 
            stdtheta= TRAp.MM_LHCP(1); bias=TRAp.MM_LHCP(2) ; %%%theta in degrees%%%
            prL(:,:,1)=funcpol(thetaRx_pier,stdtheta^2,bias);  %exp MM Left pol antenna
            prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);
        end
        
    case 'deg'
        
        if flagTxMismatch == 0
            pt(1) = 1;                                           %Transmitted ideal Right pol.vector
            pt(2) = complex(0,-1);
        else
            pt(1) = cosd(thetaTx_pier);                           %Right pol.vector with mismatch
            pt(2) = complex(0,-1);
        end
        
        if flagRxmismatch == 0
            prR(:,:,1) = ones(Ny_surf, Nx_surf);                     %Received ideal Right pol.vector
            prR(:,:,2) = complex(0,-1).*ones(Ny_surf, Nx_surf);
            prL(:,:,1) = ones(Ny_surf, Nx_surf);                     %Received ideal Left pol.vector
            prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);
        else
%             prR(:,:,1) = cosd(thetaRx_pier);                          %Right pol.vector with mismatch
%             prR(:,:,2) = complex(0,-1).*ones(Ny_surf, Nx_surf);
%             prL(:,:,1) = cosd(thetaRx_pier);                          %Left pol.vector with mismatch
%             prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);            
            stdtheta= TRAp.MM_RHCP(1); bias=TRAp.MM_RHCP(2) ; %%%theta in degrees%%%
            prR(:,:,1)=funcpol(thetaRx_pier,stdtheta^2,bias);  %exp MM Right pol antenna
            prR(:,:,2)=complex(0,-1).*ones(Ny_surf, Nx_surf); 
            stdtheta= TRAp.MM_LHCP(1); bias=TRAp.MM_LHCP(2) ; %%%theta in degrees%%%
            prL(:,:,1)=funcpol(thetaRx_pier,stdtheta^2,bias);  %exp MM Left pol antenna
            prL(:,:,2) = complex(0,1).*ones(Ny_surf, Nx_surf);

        end
end