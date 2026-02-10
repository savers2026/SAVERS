%%*****************************************************************
%   NTCMproj File: runNTCM.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%   This module output the vTEC(s), sTEC(s), and Ionospheric delay(s) computed for each row of the input data matrix.
%   It calls all the auxiliary functions and the NTCM G implementation.
%
%   CONSTANTS:
%   --------
%   pi_gal = 3.1415926535898  => Ratio of a circle's circumference to its diameter (see Table 2)
%
%   INPUT:
%   --------
%   inputData => Matrix containign for each row the following paramters: 
%     - Column Index:  [ 1  | 2   | 3   | 4   | 5   | 6            | 7     | 8         | 9            | 10          | 11             ]
%     - Column Param   [ai0 | ai1 | ai2 | DoY | UTC | Rx-longitude | Rx-latitude | Rx-Height | SV-longitude | SV-latitude | SV-Height] 
%     where:
%      - ai0 is the Effective Ionisation Level 1st order parameter [sfu]
%      - ai1 is the Effective Ionisation Level 2nd order parameter [sfu/deg]
%      - ai2 is the Effective Ionisation Level 3rd order parameter[sfu/deg]
%      - DoY is the Day of Year [dimensionless]
%      - UTC is the Universal time [hours]
%      - Rx-longitude is the User receiver Geodetic latitude [deg]
%      - Rx-latitude is the User receiver Geodetic longitude [deg]
%      - Rx-Height is the User receiver Geodetic height [meters]
%      - SV-longitude is the Satellite Geodetic latitude [deg]
%      - SV-latitude is the Satellite Geodetic longitude [deg]
%      - SV-Height is the Satellite Geodetic height Range [meters]
%
%   carrFreq => Carrier Frequency [Hz]
%
%   OUTPUT:
%   --------
%   vTEC => Vertical TEC [TECU]
%   sTEC => Slant TEC [TECU]
%   IonoDelay_m => Ionospheric delay for the input carrier frequency (carrFreq) [m]
%
%   REFERENCE:
%   --------
%	[1] European GNSS (Galileo) Open Service - NTCM G Ionospheric Model Description, Issue 1.0, European Commission (EC)
%   [2] NTCM G Software Package User Guide
% ******************************************************************
%%
function [vTEC, sTEC, IonoDelay_m] = runNTCM(inputData, carrFreq)

% ---------------------
% Init Constants
% ---------------------
pi_gal = 3.1415926535898;

% ---------------------
% Check range of input parameters
% ---------------------
[isValid, isCarrValid] = checkRanges(inputData,carrFreq);

% ---------------------
% Read input parameters
% ---------------------
brdcIonoParam = inputData(:,1:3);
matrixParam   = inputData(:,4:end);

% ---------------------
% Start Procedure
% ---------------------
nrOfExc = length(matrixParam(:,1));
deg2rad = pi_gal/180;

% Init output vectors
vTEC        = nan(nrOfExc,1);
sTEC        = nan(nrOfExc,1);
IonoDelay_m = nan(nrOfExc,1);

if ~isCarrValid
    carrFreq = nan;
end


% disp(['> ',num2str(sum(isValid),'%d'),' valid set(s) found'])
% disp('> NTCM is running...')

% Start loop 
for itc = 1 : nrOfExc
    % If any of the input parameters is not within the range, skip the Iono
    % estimation for that set
    if ~isValid
        continue;
    end
    doy     = matrixParam(itc,1);
    UTC     = matrixParam(itc,2);
    % Get User position in Geodetic coordinates
    rxPos.llh_deg  = [matrixParam(itc,4),matrixParam(itc,3),matrixParam(itc,5)];
    % Get Satellite position in Geodetic coordinates
    svPos.llh_deg  = [matrixParam(itc,7),matrixParam(itc,6),matrixParam(itc,8)];
    rxPos.llh_rad  = [rxPos.llh_deg(1).*deg2rad,rxPos.llh_deg(2).*deg2rad,rxPos.llh_deg(3)];
    svPos.llh_rad  = [svPos.llh_deg(1).*deg2rad,svPos.llh_deg(2).*deg2rad,svPos.llh_deg(3)];
    
    % *********************************************
    % STEP 1: convert Geodetic to Cartesian
    svPos.xyz  = llh2xyz(svPos.llh_rad);
    rxPos.xyz  = llh2xyz(rxPos.llh_rad);
    % *********************************************
    % STEP 2: get Azimuth and Elevation
    relVec     = svPos.xyz - rxPos.xyz; % ******   Eq. 21
    iDoa       = computeDoA(rxPos.llh_rad,relVec);
    % *********************************************
    % STEP 3: get Iono Pierce Point (IPP)
    [ippCoord]      = getIonoPiercePoint(rxPos.llh_rad,iDoa.ele_rad,iDoa.azi_rad);
    % *********************************************
    % STEP 4: convert to Local Time
    timeZone   = -((ippCoord(2)/deg2rad) /15); % ******  Eq.28 - First part
    % apply local time-zone
    LT = (UTC - timeZone);                % ******  Eq.28 - Second part
    % *********************************************
    % STEP 5: Run NTCM
    % Check elevation
    if iDoa.ele_rad<0
        vTEC(itc) = nan;
        sTEC(itc) = nan;
        continue;
    else
        % Run NTCM-G and obtain vertical TEC [TECU]
        [vTEC(itc)]      = NTCM_G(brdcIonoParam(itc,:), doy, LT, ippCoord);
        % Obtain the Mapping Function (MF)  
        % Re=6371Km is the Earth mean radius, and Hi=450Km is the IPP height
        MF               = 1./ sqrt(1-(0.934027268728925.*sin(0.9782*(pi_gal/2-iDoa.ele_rad))).^2);       % ******  Eq.33
        % Obtain slant TEC [TECU]
        sTEC(itc)        = vTEC(itc) .* MF;                                % ******  Eq.16
        % Convert slant TEC to propagation delay expressed in [m]
        IonoDelay_m(itc) = sTEC(itc) .*(40.3.*1e16)./carrFreq^2;           % ******  Eq.1
    end
end

%%*****************************************************************
%   NTCMproj File: checkRanges.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%   This module checks if the input parameters are within the ranges as
%   specified in [2] NTCM G Software Package User Guide
%
%   CONSTANTS:
%   --------
%   ai0_range => Effective Ionisation Level 1st order parameter, Range	[0, 512]sfu
%   ai1_range => Effective Ionisation Level 2nd order parameter	Range [-4, 4]sfu/deg
%   ai2_range => Effective Ionisation Level 3rd order parameter	Range [-0.25, 0.25]sfu/deg
%   doy_range => Day of Year Range [1, 366] dimensionless
%   utc_range => Universal time [0,24]hours
%   rlt_range => User receiver Geodetic latitude Range [-90, 90]deg
%   rln_range => User receiver Geodetic longitude Range [-180, 180]deg
%   rht_range => User receiver Geodetic height Range [-4000, 400000]meters
%   slt_range => Satellite Geodetic latitude Range [-90, 90]deg
%   sln_range => Satellite Geodetic longitude Range [-180, 180]deg
%   sht_range => Satellite Geodetic height Range [450000, 60000000]meters
%   crr_range => Carrier Frequency Range [0, 20000e6]Hz
%
%   INPUT:
%   --------
%   inputData => Matrix containign for each row the following paramters: 
%     - Column Index:  [ 1  | 2   | 3   | 4   | 5   | 6            | 7     | 8         | 9            | 10          | 11             ]
%     - Column Param   [ai0 | ai1 | ai2 | DoY | UTC | Rx-longitude | Rx-latitude | Rx-Height | SV-longitude | SV-latitude | SV-Height] 
%   carrFreq => Signal Carrier Frequency [Hz]
%
%   OUTPUT:
%   --------
%   isValid => Vector indicating for each row of <inputData> if the parameters are within the ranges
%   isCarrValid => Flag indicating if the Carrier frequency is within the range
%
%   REFERENCE:
%   --------
%   [2] NTCM G Software Package User Guide
% ******************************************************************

function [isValid, isCarrValid] = checkRanges(inputData, carrFreq)

ai0_range = [0,512];       
ai1_range = [-4,4];
ai2_range = [-0.25,0.25];
doy_range = [1,366];
utc_range = [0,24];
rlt_range = [-90,90];
rln_range = [-180,180];
rht_range = [-4e3,400e3];
slt_range = [-90,90];
sln_range = [-180,180];
sht_range = [450e3,60000e3];
crr_range = [0,20000e6];

% Check input ranges
isValid = (inputData(:,1) >=ai0_range(1) & inputData(:,1) <=ai0_range(2)) & ...
          (inputData(:,2) >=ai1_range(1) & inputData(:,2) <=ai1_range(2)) & ...
          (inputData(:,3) >=ai2_range(1) & inputData(:,3) <=ai2_range(2)) & ...
          (inputData(:,4) >=doy_range(1) & inputData(:,4) <=doy_range(2)) & ...
          (inputData(:,5) >=utc_range(1) & inputData(:,5) <=utc_range(2)) & ...
          (inputData(:,6) >=rln_range(1) & inputData(:,6) <=rln_range(2)) & ...
          (inputData(:,7) >=rlt_range(1) & inputData(:,7) <=rlt_range(2)) & ...
          (inputData(:,8) >=rht_range(1) & inputData(:,8) <=rht_range(2)) & ...
          (inputData(:,9) >=sln_range(1) & inputData(:,9) <=sln_range(2)) & ...
          (inputData(:,10)>=slt_range(1) & inputData(:,10)<=slt_range(2)) & ...
          (inputData(:,11)>=sht_range(1) & inputData(:,11)<=sht_range(2));
      
isCarrValid = (carrFreq >=crr_range(1) & carrFreq <=crr_range(2));
end

%%*****************************************************************
%   NTCMproj File: computeDoA.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%    This module computes the satellite azimuth and elevation angles (Described in Section 2.5.2)    
%
%   CONSTANTS:
%   --------
%   pi_gal = 3.1415926535898  => Ratio of a circle's circumference to its diameter (see Table 2)
%
%   INPUT:
%   --------
%   llhUserRad => Receiver position in  Latitude[rad], Longitude[rad], height[m]
%   losVec     => Line of Sight (LoS) vector
%
%   OUTPUT:
%   --------
%   svDoa.ele_rad      => Elevation (DoA) [rad]
%   svDoa.azi_rad      => Azimuth (DoA)   [rad]
%
%   REFERENCE:
%   --------
%	[1] European GNSS (Galileo) Open Service - NTCM G Ionospheric Model Description, Issue 1.0, European Commission (EC)
% ******************************************************************
%%
function svDoa = computeDoA(llhUserRad, losVec)
% 
pi_gal = 3.1415926535898; % Ratio of a circle's circumference to its diameter (see Table 2)

phi = llhUserRad(:,1);
lam = llhUserRad(:,2);
sinphi = sin(phi);
cosphi = cos(phi);
sinlam = sin(lam);
coslam = cos(lam);
cosphicoslam = cosphi .* coslam;
cosphisinlam = cosphi .* sinlam;
sinphicoslam = sinphi .* coslam;
sinphisinlam = sinphi .* sinlam; 

% Compute rotation matrix (Eq. 21)
dV_rot = [-sinphicoslam -sinphisinlam cosphi ;...
          -sinlam        coslam       0;...
           cosphicoslam  cosphisinlam sinphi] * losVec';
       
% Calculate azimuth of DoA (Eq. 22)
svDoa.azi_rad  = atan2(dV_rot(2),dV_rot(1));
if svDoa.azi_rad<0
    svDoa.azi_rad = svDoa.azi_rad + 2*pi_gal;
end       

% Calculate elevation of DoA (Eq. 23)      
svDoa.ele_rad  = ( 0.5 - atan2( sqrt(dV_rot(1).^2+dV_rot(2).^2), dV_rot(3)) ./ pi_gal ) .* pi_gal;

end

%%*****************************************************************
%   NTCMproj File: getIonoPiercePoint.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%   This module computes the geographic latitude and longitude of the Ionospheric Pierce Point (IPP) (Described in Section 2.5.3). 
%   IPP is the point where the Line of Sight (LoS) intersects the reference ionospheric layer, this latter set at an altitude of 450Km.
%
%   CONSTANTS:
%   --------
%   Re     = 6371.0*1e3       => Earth mean radius (see Table 2)
%   h_IPP  = 450e3            => Ionospheric Pierce Point height (see Table 2)
%   pi_gal = 3.1415926535898  => Ratio of a circle's circumference to its diameter (see Table 2)
%
%   INPUT:
%   --------
%   llhUser    => Receiver position in  Latitude[rad], Longitude[rad], height[m]
%   ElRad      => Elevation (DoA) [rad]
%   AzRad      => Azimuth (DoA)   [rad]
%
%   OUTPUT:
%   --------
%   IppCoord => Pierce Point coordinate in  Latitude[rad], Longitude[rad]   
%
%
%   REFERENCE:
%   --------
%	[1] European GNSS (Galileo) Open Service - NTCM G Ionospheric Model Description, Issue 1.0, European Commission (EC)
% ******************************************************************
%% 
function [IppCoord] = getIonoPiercePoint(llhUser,ElRad,AzRad)
Re        = 6371.0*1e3;        % Earth mean radius (see Table 2)
h_IPP     = 450e3;             % Ionospheric Pierce Point height (see Table 2)
pi_gal    = 3.1415926535898;   % Ratio of a circle's circumference to its diameter (see Table 2)

r = h_IPP + Re;

nSmp      = length(ElRad);
IppCoord  = nan(length(ElRad),2);

for ipos = 1 : nSmp
    % Calculate Earth's central angle (Eq. 24)
    Psi_pp            = pi_gal/2 - ElRad(ipos) - asin(Re/r.*cos(ElRad(ipos)));
    sinPsi_pp         = sin(Psi_pp);
    % Calculate pierce Point latitude (Eq. 25)
    IppCoord(ipos,1)  = asin(sin(llhUser(1)).*cos(Psi_pp)+cos(llhUser(1)).*sinPsi_pp.*cos(AzRad(ipos)));
    % Calculate pierce Point longitude (Eq. 26)
    IppCoord(ipos,2)  = llhUser(2) + asin(sinPsi_pp.*sin(AzRad(ipos))./cos(IppCoord(ipos,1)));       
end
end
%%*****************************************************************
%   NTCMproj File: llh2xyz.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%   This module converts the geodetic WGS84 coordinates to ECEF (Described in Sec.2.5.1)
%
%   CONSTANTS:
%   --------
%   a  = 6378137.0         => Semi-major axis [m]
%   b  = 6356752.3142      => Semi-major axis [m]
%   e2 = 0.006694380004261 => Eccentricity (e) of the ellipsoid squared (e2=1-b^2/a^2)
%
%   INPUT:
%   --------
%   llh => WGS84 coordinates in Latitude[rad], Longitude[rad], height[m]
%
%   OUTPUT:
%   --------
%   xyz => ECEF coordinates in x[m], y[m], z[m]
%
%   REFERENCE:
%   --------
%	[1] European GNSS (Galileo) Open Service - NTCM G Ionospheric Model Description, Issue 1.0, European Commission (EC)
% ******************************************************************
%%
function [xyz] = llh2xyz(llh)

% Semi-major axis [m]
a  = 6378137.0;
% Eccentricity (e) of the ellipsoid squared (e2)
% (It can be also obtained as <e2 = f*(2-f)> where f is the flattening)
e2 = 0.006694380004261;

% (Eq. 19)
v = a ./ sqrt(1-e2*sin(llh(:,1)).^2);
% Compute x component (Eq. 16)
xyz(:,1) = (v+llh(:,3)) .* cos(llh(:,1)) .* cos(llh(:,2));
% Compute y component (Eq. 17)
xyz(:,2) = (v+llh(:,3)) .* cos(llh(:,1)) .* sin(llh(:,2));
% Compute z component (Eq. 18)
xyz(:,3) = (v*(1-e2)+llh(:,3)) .* sin(llh(:,1));
end

%%*****************************************************************
%   NTCMproj File: NTCM_G.m
%%*****************************************************************
%   @author      Matteo Sgammini
%   @reviewer    Francesco Menzione
%   @ingroup     NTCM_JRC
%   @copyright   Joint Research Centre (JRC), 2022
%   This software has been released as free and open source software
%   under the terms of the European Union Public Licence (EUPL), version 1
%   Questions? Submit your query at https://www.gsc-europa.eu/contact-us/helpdesk
%%*****************************************************************
%   Code generated for Matlab model 'NTCM_G'
%   Model version                  : 1.0
%   MatLab version                 : 9.7.0.1190202 (R2019b)
%
%%*****************************************************************
%   FUNCTION:
%   --------
%   This module implements the NTCM G model.
%   Output of the module is the integrated VTEC. 
%
%   CONSTANTS:
%   --------
%
%   k1,k2,..,k12 = (see Table 3) => NTCM model coefficients (Described in Section 2.4.3)
%   Pf1          = 0.4 => Phase shift for the solar zenith angle dependency [dimensionless]
%   LTd          = 14  => Time phase shift for the Local Time [hours]
%   doy_DAV      = 18  => Time phase shift for the Annual variation [days]
%   doy_DSAV     = 6   => Time phase shift for the Semi-annual variation [days]
%   phi_c1       = 16  => Northward ionisation crest [deg]
%   phi_c2       = -10 => Southward ionisation crest [deg]
%   sigma_c1     = 12  => Best-fit value for the 1st component of the Ionisation crests  [deg]
%   sigma_c2     = 13  => Best-fit value for the 2nd component of the Ionisation crests  [deg]
%   phiGNP_deg   =  79.74 => Latitude of the geomagnetic North pole [deg]
%   lamGNP_deg   = -71.78 => Longitude of the geomagnetic North pole [deg]
%   pi_gal       = 3.1415926535898  => Ratio of a circle's circumference to its diameter (see Table 2)
%
%   INPUT:
%   --------
%   brdcIonoParam  => ai0,ai1,ai2 Effective Ionisation Level parameters broadcast in the Galileo Navigation Message
%   doy            => Day of the Year
%   LT             => Local Time [hours]
%   IppCoord       => Pierce Point Coordinates in  Latitude[rad], Longitude[rad] 
%
%   OUTPUT:
%   --------
%   vTEC => Vertical TEC in [TECU]
%
%   REFERENCE:
%   --------
%	[1] European GNSS (Galileo) Open Service - NTCM G Ionospheric Model Description, Issue 1.0, European Commission (EC)
%   [3] Hoque MM, N Jakowski, R Orus Perez, "Fast ionospheric correction using Galileo Az coefficients and the NTCM model,
%       GPS Solutions, 2019, doi: 10.1007/s10291-019-0833-3. https://doi.org/10.1007/s10291-019-0833-3
% ******************************************************************
%%
function [vTEC] = NTCM_G(brdcIonoParam,doy,LT,IppCoord)

k1 = 0.92519; k2 = 0.16951; k3 = 0.00443; k4 = 0.06626; k5 = 0.00899; k6 = 0.21289; k7 = -0.15414; k8 = -0.38439; k9 = 1.14023; k10 = 1.20556; k11 = 1.41808; k12 = 0.13985; % rms
Pf1        = 0.4;
LTd        = 14;
doy_DAV    = 18;
doy_DSAV   = 6;
phi_c1     = 16;     %degrees
phi_c2     = -10;    %degrees
sigma_c1   = 12;     %degrees
sigma_c2   = 13;     %degrees
phiGNP_deg =  79.74; %degrees
lamGNP_deg = -71.78; %degrees
pi_gal     = 3.1415926535898;   % Ratio of a circle's circumference to its diameter (see Table 2)

deg2rad    = pi_gal/180;

%Geomagnetic North Pole coordinate
phiGNP     = phiGNP_deg*deg2rad; % latitude in radian
lamGNP     = lamGNP_deg*deg2rad; % longitude in radian 
% Proxy measure of the solar activity level
azPar      = sqrt(          brdcIonoParam(1).^2 + ...
                  1633.33 * brdcIonoParam(2).^2 + ...
                  4802000 * brdcIonoParam(3).^2 + ...
                  3266.67 * brdcIonoParam(1) .* brdcIonoParam(3) );                              % ******  Eq.2

phiIPP = IppCoord(:,1);
lamIPP = IppCoord(:,2);

% Sun's declination [rad]  
delta = 23.44*sin(0.9856*(doy-80.7)*deg2rad)*deg2rad;                                            % ******  Eq.28

cos_lat_delta    = cos(phiIPP-delta);
solzenith_factor = cos_lat_delta + Pf1;                                                          % ******  Eq.29
cos_solzenith    = cos_lat_delta - phiIPP/(pi_gal/2).*sin(delta);                                % ******  Eq.30

latm_rad = asin(sin(phiIPP).*sin(phiGNP)+cos(phiIPP).*cos(phiGNP).*cos(lamIPP-lamGNP));          % ******  Eq.31
latm_deg = latm_rad/deg2rad; 

Vd     = 2*pi_gal*(LT-LTd)/24;                                                                   % ******  Eq.5
Vsd    = 2*pi_gal*LT/12;                                                                         % ******  Eq.6
Vtd    = 2*pi_gal*LT/8;                                                                          % ******  Eq.7
Va     = 2*pi_gal*(doy-doy_DAV)/365.25;                                                          % ******  Eq.9
Vsa    = 4*pi_gal*(doy-doy_DSAV)/365.25;                                                         % ******  Eq.10

cosVd  = cos(Vd);
cosVsd = cos(Vsd);
sinVsd = sin(Vsd);
cosVtd = cos(Vtd);
sinVtd = sin(Vtd);

EC1 = -(latm_deg-phi_c1)^2/2/sigma_c1^2;                                                         % ******  Eq.13
EC2 = -(latm_deg-phi_c2)^2/2/sigma_c2^2;                                                         % ******  Eq.14
exp_EC1 = exp(EC1);
exp_EC2 = exp(EC2);

% F1: Local time dependency
F1 = solzenith_factor + cos_solzenith.*(k1*cosVd+k2*cosVsd+k3*sinVsd+k4*cosVtd+k5*sinVtd);       % ******  Eq.4
% F2: Seasonal Dependency
F2 = 1 + k6.*cos(Va)+k7.*cos(Vsa);                                                               % ******  Eq.8;
% F3: Geomagnetic field dependency
F3 = 1 + k8.*cos(latm_rad);                                                                      % ******  Eq.11
% F4: Equatorial anomaly dependency
F4 = 1 + k9.*exp_EC1+k10.*exp_EC2;                                                               % ******  Eq.12
% F5: Solar activity dependency
F5 = k11 + k12.*azPar;                                                                           % ******  Eq.15

% Compute the vertical TEC
vTEC = F1 .* F2 .* F3 .* F4 .* F5;                                                               % ******  Eq.3
end

% disp('Finished!')
end