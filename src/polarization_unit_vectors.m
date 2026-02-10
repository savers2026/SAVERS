%% ************************************************************************
% MODULE NAME:      polarization_unit_vectors.m
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

function [pt, prR, prL] = polarization_unit_vectors(flagTxMismatch, flagRxmismatch, thetaTx_pier, thetaRx_pier, tag)


switch tag
    case 'rad'
        
        if flagTxMismatch == 0
            pt(1) = 1;                   %Transmitted ideal Right pol.vector
            pt(2) = complex(0,-1);
        else
            pt(1) = cos(thetaTx_pier);   %Right pol.vector with mismatch
            pt(2) = complex(0,-1);
        end
        
        if flagRxmismatch == 0
            prR(1) = 1;                  %Received ideal Right pol.vector
            prR(2) = complex(0,-1);
            prL(1) = 1;                  %Received ideal Left pol.vector
            prL(2) = complex(0,1);
        else
            prR(1) = cos(thetaRx_pier);  %Right pol.vector with mismatch
            prR(2) = complex(0,-1);
            prL(1) = cos(thetaRx_pier);  %Left pol.vector with mismatch
            prL(2) = complex(0,1);
            
        end
        
    case 'deg'
        
        if flagTxMismatch == 0
            pt(1) = 1;                   %Transmitted ideal Right pol.vector
            pt(2) = complex(0,-1);
        else
            pt(1) = cosd(thetaTx_pier); %Right pol.vector with mismatch
            pt(2) = complex(0,-1);
        end
        
        if flagRxmismatch == 0
            prR(1) = 1;                  %Received ideal Right pol.vector
            prR(2) = complex(0,-1);
            prL(1) = 1;                  %Received ideal Left pol.vector
            prL(2) = complex(0,1);
        else
            prR(1) = cosd(thetaRx_pier); %Right pol.vector with mismatch
            prR(2) = complex(0,-1);
            prL(1) = cosd(thetaRx_pier); %Left pol.vector with mismatch
            prL(2) = complex(0,1);
            
        end
        
        
        
        
end