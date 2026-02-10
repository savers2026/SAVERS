%% ************************************************************************
% MODULE NAME:      transformation_abs_flr_mat.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Guerriero and Davide Comite
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      Sep 2019 by L. Guerriero
% Main updates:          2021 by Davide Comite
% Released to ESA:       December 2022
%% ************************************************************************
% FUNCTION
% The module converts the incidence and scattering angles from  
% the Global SAVERS reference system to the facet local reference (FLR) at  
% each pixel point (xp,yp,zp), that is the ENU rotated frame which has been 
% also tilted to take into account of the slope and aspect of the facet.
%% ************************************************************************
% REFERENCE
% HydroGNSS E2E Simulator User Guide 
% HydroGNSS E2E Simulator Algorithm Theoretical Baseline Document
%% ************************************************************************

function [thetaprimeloc,phiprimeloc,thetasprimeloc,phisprimeloc] = transformation_abs_flr_mat(beta,alpha,thetai,thetas,phis)

Ny_surf = size(thetas,1);
Nx_surf = size(thetas,2);

U11 =  cosd(beta).*cosd(alpha);  %questa rotazione di beta dà rotazione in verso antiorario nel piano zx (va verso x<0)
U12 =  sind(alpha).*cosd(beta);
U13 = -sind(beta);
U21 = -sind(alpha);
U22 =  cosd(alpha);
U23 =  0;
U31 =  cosd(alpha).*sind(beta);
U32 =  sind(beta).*sind(alpha);
U33 =  cosd(beta);

%Principale

incprimeloc(:,:,1) = U11.*sind(thetai) + U12.*0 - U13.*cosd(thetai); 
incprimeloc(:,:,2) = U21.*sind(thetai) + U22.*0 - U23.*cosd(thetai); 
incprimeloc(:,:,3) = U31.*sind(thetai) + U32.*0 - U33.*cosd(thetai); 

zprimeloc(:,:,1) = zeros(Ny_surf,Nx_surf);
zprimeloc(:,:,2) = zeros(Ny_surf,Nx_surf);
zprimeloc(:,:,3) = ones(Ny_surf,Nx_surf);

dotzetainc   = dot(zprimeloc,incprimeloc,3);

% if (1-abs(dotzetainc))<1.0e-10 
%     if dotzetainc < 0.0 dotzetainc = -1.0; end
%     if dotzetainc > 0.0 dotzetainc = 1.0; end
% end

thetaprimeloc = real(acosd(dotzetainc));
thetaprimeloc = 180.0-thetaprimeloc;

xprimeloc(:,:,1) = ones(Ny_surf,Nx_surf);
xprimeloc(:,:,2) = zeros(Ny_surf,Nx_surf);
xprimeloc(:,:,3) = zeros(Ny_surf,Nx_surf);

cosphiprimeloc = dot(xprimeloc,incprimeloc,3)./sind(thetaprimeloc);

% if (1.0-abs(cosphiprimeloc))<1.0e-10
%     if cosphiprimeloc<0.0 cosphiprimeloc = -1.0; end
%     if cosphiprimeloc>0.0 cosphiprimeloc = 1.0; end
% end

yprimeloc(:,:,1) = zeros(Ny_surf,Nx_surf);
yprimeloc(:,:,2) = ones(Ny_surf,Nx_surf);
yprimeloc(:,:,3) = zeros(Ny_surf,Nx_surf);

phiprimeloc  = real(acosd(cosphiprimeloc));
phiprimelocb = real(asind(dot(yprimeloc,incprimeloc,3)./sind(thetaprimeloc)));

% if phiprimelocb<0 
%     phiprimeloc = -phiprimeloc; 
% end

ind000 = phiprimelocb < 0;
phiprimeloc(ind000) = -phiprimeloc(ind000);

% if sind(thetaprimeloc) == 0   
%     phiprimeloc = 0; 
% end

ind001 = sind(thetaprimeloc) == 0   ;
phiprimeloc(ind001) = 0;

% sprimeloc=U*s;       %vett. s espresso nel sist. locale

sprimeloc(:,:,1) = U11.*sind(thetas).*cosd(phis) + U12.*sind(thetas).*sind(phis) + U13.*cosd(thetas); 
sprimeloc(:,:,2) = U21.*sind(thetas).*cosd(phis) + U22.*sind(thetas).*sind(phis) + U23.*cosd(thetas); 
sprimeloc(:,:,3) = U31.*sind(thetas).*cosd(phis) + U32.*sind(thetas).*sind(phis) + U33.*cosd(thetas);

thetasprimeloc = real(acosd(dot(zprimeloc,sprimeloc,3)));

phisprimeloc  = real(acosd(dot(xprimeloc,sprimeloc,3)./sind(thetasprimeloc)));
phisprimelocb = real(asind(dot(yprimeloc,sprimeloc,3)./sind(thetasprimeloc)));

% if phisprimelocb<0 
%     phisprimeloc=-phisprimeloc; 
% end
% 
ind002 = phisprimelocb < 0 ;
phisprimeloc(ind002) = -phisprimeloc(ind002);
% 
% if sind(thetasprimeloc)==0.0 
%     phisprimeloc=0.0; 
% end
% 
ind003 = sind(thetasprimeloc) == 0;
phisprimeloc(ind003) = -phisprimeloc(ind003);
% 
% end
