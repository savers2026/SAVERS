%--------------------------------------------------------------------------
% This Function transforms a vector from ENU rotated (whih is, the
% SAVERS Global Reference system) to EFEC
%------------------------------------------------------------------------

% Written ex-novo: Nov 2021 for HSAVERS Davide Comite


function pECEF = local2ecef_dc_mat(phi_ENU2local, pLOCAL, OenuINecef) 


pENU  = rotation_counterclockwise_mat(pLOCAL, phi_ENU2local, 'rad');

pECEF = enu2ecef_dc_mat(pENU, OenuINecef);  


end
