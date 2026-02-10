% ECEF2ENU conversion

% Convert from Earth Centered Eart Fixed to local East, North, Up (ENU) coordinates
% Notice: since Earth is an ellispoid, East is not exactly oriented like longitude,
% wheras north is vs latitude

%%%%Input
%   p:           Vector to be converted in ENU - column vector
%   OENUinECEF:  Origin of ENU in ECEF, i.e., the SP in ECEF
%   type:        A string, v, if present indicates conversion of velocity vector

%%%%Output
%   Penu: colum vectors in ENU coordinates

% Written ex-novo: Oct 2021 for HSAVERS Davide Comite


function [Penu] = ecef2enu_dc_mat(vectECEF, OENUinECEF, type) 

%% Input Checks

%check if vect is column vector
%if isrow(vectECEF), vectECEF = vectECEF.'; end

%% Get lat and lon for reference O

[ph, lm, hs] = xyz2llh(OENUinECEF(1), OENUinECEF(2), OENUinECEF(3));

%% Compute rotation matrix: from ecef 2 enu

%B_ECEF2ENU = [ -sin(lm)               cos(lm)       0   ; 
%               -sin(ph)*cos(lm)  -sin(ph)*sin(lm) cos(ph);...
%                cos(ph)*cos(lm)   cos(ph)*sin(lm) sin(ph)];

  
   
%% Conversion to ENU is done first by translation to Oecef and then rotation

%Translation - not applicable for velocity vector
if nargin == 2

          %pt = vectECEF - OENUinECEF;
          pt(:,:,1) = vectECEF(:,:,1) - OENUinECEF(1);
          pt(:,:,2) = vectECEF(:,:,2) - OENUinECEF(2);
          pt(:,:,3) = vectECEF(:,:,3) - OENUinECEF(3);

          %Rotation
          %Penu = B_ECEF2ENU*pt;

          Penu(:,:,1) =  -pt(:,:,1).*sin(lm)         + pt(:,:,2).*cos(lm)         +        0;
          Penu(:,:,2) =  -pt(:,:,1).*sin(ph)*cos(lm) - pt(:,:,2).*sin(ph)*sin(lm) +  pt(:,:,3).*cos(ph);
          Penu(:,:,3) =   pt(:,:,1).*cos(ph)*cos(lm) + pt(:,:,2).*cos(ph)*sin(lm) +  pt(:,:,3).*sin(ph);

elseif nargin == 3

          %Penu = B_ECEF2ENU*vectECEF;

          Penu(:,:,1) =  -pt(:,:,1).*sin(lm)         + pt(:,:,2).*cos(lm)         +        0;
          Penu(:,:,2) =  -pt(:,:,1).*sin(ph)*cos(lm) - pt(:,:,2).*sin(ph)*sin(lm) +  pt(:,:,3).*cos(ph);
          Penu(:,:,3) =   pt(:,:,1).*cos(ph)*cos(lm) + pt(:,:,2).*cos(ph)*sin(lm) +  pt(:,:,3).*sin(ph);

end



