function [phi, lambda, h] = xyz2llh_mat(pLocal)

%Convert ECEF point into lat (rad), lon (rad), h (m)

% Written ex-novo: Feb 2022 for HSAVERS Davide Comite

%Based on: J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates ...
%           to geodetic coordinates," Aerospace and Electronic Systems, ...
%           IEEE Transactions on, vol. 30, pp. 957-961, 1994.

  %Earth semimajor axis in meters
  a = 6378137; 
  %Reciprocal flattening
  f = 1/298.257223563; 
  %Semi-minor axis
  b = a*(1-f);
 
  %First eccentricity squared
  e2  = 2*f - f^2;
  %Second eccentricity squared
  ep2 = f*(2-f)/((1-f)^2); 
 
  r2 = pLocal(:,:,1).^2 + pLocal(:,:,2).^2;
  r  = sqrt(r2);
  
  E2 = a^2 - b^2;
  F  = 54.*b^2.*pLocal(:,:,3).^2;

  G  = r2 + (1-e2).*pLocal(:,:,3).^2 - e2.*E2;
  c  = (e2.*e2.*F.*r2)./(G.*G.*G);
  s  = ( 1 + c + sqrt(c.*c + 2*c) ).^(1/3);
  P  = F./(3*(s+1./s+1).^2.*G.*G);
  Q  = sqrt(1+2.*e2.*e2.*P);
  
  ro  = -(e2*P.*r)./(1 + Q) + sqrt((a.*a/2)*(1 + 1./Q) - ((1-e2)*P.*pLocal(:,:,3).^2)./(Q.*(1 + Q)) - P.*r2./2);
  tmp = (r - e2.*ro).^2;
  
  U   = sqrt(tmp + pLocal(:,:,3).^2);
  V   = sqrt(tmp + (1-e2).*pLocal(:,:,3).^2);

  zo  = (b.^2.*pLocal(:,:,3))./(a.*V);
 
  h      = U.*(1 - b.^2./(a.*V));
  phi    = atan((pLocal(:,:,3) + ep2.*zo)./r);
  lambda = atan2(pLocal(:,:,2), pLocal(:,:,1));
