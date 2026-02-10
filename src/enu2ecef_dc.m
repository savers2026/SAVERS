% ENU2ECEF conversion

% Convert to Earth Centered Eart Fixed
% from local East, North, Up (ENU) coordinates
%
% Input:
%   V_enu:  [3 N] matrix with column vectors to be converted, N >= 1
%   O_ecef:  [3 1] column vector with ecef coordinates of ref point O
%           (the one used for origin in the ENU reference)
% Output:
%   P_ecef: [3 N] matrix with colum vectors in ECEF coordinates

% Written ex-novo: Oct 2021 for HSAVERS Davide Comite


function [P_ecef] = enu2ecef_dc(vectENU, O_ecef) 

%% Input Checks

%check if vect is column vector
if isrow(vectENU), vectENU = vectENU.'; end

%% Get lat and lon for reference O 

[ph, lm, hs] = xyz2llh(O_ecef(1), O_ecef(2), O_ecef(3) );

%% Compute rotation matrix: from ecef 2 enu

Bmat_ECEF2ENU = [ -sin(lm)               cos(lm)       0; 
                  -sin(ph)*cos(lm) -sin(ph)*sin(lm) cos(ph);...
                   cos(ph)*cos(lm)  cos(ph)*sin(lm) sin(ph)];
               
Bmat_ENU2ECEF = Bmat_ECEF2ENU.';               
  
%% Conversion to ENU is done first by roation into ECEF then Translation
% Rotation
Vecef = Bmat_ENU2ECEF*vectENU;

% Translation
P_ecef = Vecef + O_ecef;

