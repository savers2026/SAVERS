%% ************************************************************************
% MODULE NAME:      transformation_abs_flr.m
% SOFTWARE NAME:    HSAVERS
% SOFTWARE VERSION: 2.5
%% ************************************************************************
% AUTHORS:      L. Guerriero
% @copyright:   Tor Vergata University of Rome and Sapienza University of Rome
% Original version:      Sep 2019 by L. Guerriero
% Main updates:          
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

function [thetaprimeloc,phiprimeloc,thetasprimeloc,phisprimeloc] = transformation_abs_flr(beta,alpha,thetai,thetas,phis)

%Vettori colonna: prod. scalare= vet'*vet, oppure dot(vet,vet)
%prodotto vettoriale: cross(vet1,vet2)
%Il prodotto vet*vet' dà una matrice
%Matrice M(3,3): prodotto riga x colonna=M*vet dà vettore colonna

U(1,1)=cosd(beta)*cosd(alpha);  %questa rotazione di beta dà rotazione in verso antiorario nel piano zx (va verso x<0)
U(1,2)=sind(alpha)*cosd(beta);
U(1,3)=-sind(beta);
U(2,1)=-sind(alpha);
U(2,2)=cosd(alpha);
U(2,3)=0;
U(3,1)=cosd(alpha)*sind(beta);
U(3,2)=sind(beta)*sind(alpha);
U(3,3)=cosd(beta);

%Principale
inc(1,1)=sind(thetai);
inc(2,1)=0;
inc(3,1)=-cosd(thetai);
zeta=[0; 0; 1];     %asse z sist principale
xics=[1; 0; 0];     %asse x sist principale

s(1,1)=sind(thetas)*cosd(phis);
s(2,1)=sind(thetas)*sind(phis);
s(3,1)=cosd(thetas);

%Locale
zetaprimeloc=[0.0; 0.0; 1.0];  % asse z  sist.locale
xprimeloc=[1.0; 0.0; 0.0];  % asse x sist.locale
yprimeloc=[0.0; 1.0; 0.0];  % asse y sist.locale
incprimeloc=U*inc;       %vett. inc. espresso nel sist. locale
dotzetainc=dot(zetaprimeloc,incprimeloc);
if (1-abs(dotzetainc))<1.0e-10 
    if dotzetainc<0.0 dotzetainc=-1.0; end
    if dotzetainc>0.0 dotzetainc=1.0; end
end
thetaprimeloc=real(acosd(dotzetainc));
thetaprimeloc=180.0-thetaprimeloc;

cosphiprimeloc=dot(xprimeloc,incprimeloc)/sind(thetaprimeloc);
if (1.0-abs(cosphiprimeloc))<1.0e-10
    if cosphiprimeloc<0.0 cosphiprimeloc=-1.0; end
    if cosphiprimeloc>0.0 cosphiprimeloc=1.0; end
end
phiprimeloc=real(acosd(cosphiprimeloc));
phiprimelocb=real(asind(dot(yprimeloc,incprimeloc)/sind(thetaprimeloc)));
if phiprimelocb<0 phiprimeloc=-phiprimeloc; end
if sind(thetaprimeloc)==0.0  phiprimeloc=0; end

sprimeloc=U*s;       %vett. s espresso nel sist. locale
thetasprimeloc=real(acosd(dot(zetaprimeloc,sprimeloc)));

phisprimeloc=real(acosd(dot(xprimeloc,sprimeloc)/sind(thetasprimeloc)));
phisprimelocb=real(asind(dot(yprimeloc,sprimeloc)/sind(thetasprimeloc)));
if phisprimelocb<0 phisprimeloc=-phisprimeloc; end
if sind(thetasprimeloc)==0.0 phisprimeloc=0.0; end
return
