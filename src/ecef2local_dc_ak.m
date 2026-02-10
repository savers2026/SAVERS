function pLOCAL = ecef2local_dc_ak(phi_ENU2local, pECEF, OenuINecef, varargin)


    if nargin == 3                     % position-like: translate
        pENU = ecef2enu_dc(pECEF, OenuINecef);
    elseif nargin == 4                 % velocity-like: skip translation
        pENU = ecef2enu_dc(pECEF, OenuINecef, 'v');
    else
        error('ecef2local_dc expects 3 or 4 input arguments.');
    end

    % Final rotation about the local Up axis
    pLOCAL = rotation_clockwise(pENU, phi_ENU2local, 'rad');
end
