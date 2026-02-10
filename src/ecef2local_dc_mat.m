%--------------------------------------------------------------------------
% This Function transforms a vector from ECEF to ENU rotated (whih is, the
% SAVERS Global Reference system)
%--------------------------------------------------------------------------

% Written ex-novo: Nov 2021 for HSAVERS Davide Comite

%--------------------------------------------------------------------------

function [pLOCAL] = ecef2local_dc_mat(phi_ENU2local, pECEF, OenuINecef) 
      
    pENU     = ecef2enu_dc_mat(pECEF, OenuINecef);  
      
    pLOCAL   = rotation_clockwise_mat(pENU, phi_ENU2local, 'rad');     

end
