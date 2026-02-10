function beta_c_opt = findBetaC_New(R0, Rs, lambda, Pt, sigma_z, ...
                                 TRApGain_rx, Ny_surf, Nx_surf, TRApGain_tx, sigmaParam_Betac, TRAp_Betac, SMvol,gammaParam,surfParam)

Pt = db2pow(Pt);
sigma_z_cm = sigma_z;
sigma_z = sigma_z / 100;

theta_i = TRAp_Betac.incAngle_enuRot;
phi_i = 0;

theta_s = sigmaParam_Betac.TH_s; %matrix
phi_s = sigmaParam_Betac.PH_s; %matrix

% D_rx = @(theta_deg, phi_deg) getGainMatrix(TRApGain_rx, theta_deg, phi_deg);
% D_tx = @(theta_deg, phi_deg) getGainMatrix(fliplr(TRApGain_tx), theta_deg, phi_deg);

% At SP, theta_i = theta_s
% phi_i = assumption: 0

R0_slant = R0./cosd(theta_i);
Rs_slant = Rs./cosd(theta_i);

%% fresnel coeff from the routine given by prof.
cd=costantedielettrica(SMvol./100,TRAp_Betac.TxFrequency);
refr_index=sqrt(cd);
[rH,rV]=fresnel_Matlab(1,refr_index,theta_i);
reflectioncoeff=(rV+rH)/2;

 P_IT = Pt * gammaParam.Gain_LHCP_onsurf_lin_atSP * TRAp_Betac.TxGainAtSP ...
           * lambda^2 * abs(reflectioncoeff)^2 ...
           / ((4*pi)^2 * (R0_slant + Rs_slant)^2);


theta_s_vector = sigmaParam_Betac.TH_s(:);  %flatten matrix to vector
phi_s_vector   = sigmaParam_Betac.PH_s(:);

% P_IT_vector = repmat(P_IT, N, 1);  %scalar P_IT to vector (we have N theta_s & phi_s and one P_IT)

% rmseFun = @(beta_c) sqrt(mean(arrayfun(@(th_s, ph_s, Pit) ...
%     (radarEquation(R0, Rs, lambda, sigma_z, beta_c, th_s, ph_s, Pt, D_rx, D_tx, ...
%                    Ny_surf, Nx_surf, ...
%                    TRAp_Betac.TxFrequency, SMvol, sigma_z_cm, ...
%                    R0_slant, Rs_slant, ...
%                    beta_c, beta_c, theta_i, phi_i, theta_s_vector , phi_s_vector) - Pit).^2, ...
%     theta_s_vector, phi_s_vector, P_IT_vector)));


beta_c = 0.00001:0.1:100; % beta_c is beta_2

beta_0 = 0.03; %fix the beta_0

pr = zeros(size(beta_c));
AbsE = zeros(size(beta_c));
coeff_matrix_all = cell(length(beta_c), 1);

for i = 1:length(beta_c)
    [pr_scalar, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, beta_c(i), Pt, ...
                              Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
                              SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
                              beta_0, beta_c(i), theta_i, phi_i, theta_s_vector, phi_s_vector, surfParam,TRAp_Betac,sigmaParam_Betac);
 
% for i = 1:length(beta_c)
%     [pr_scalar, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, beta_c(i), Pt, ...
%                               Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
%                               SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
%                               beta_c(i), beta_c(i), theta_i, phi_i, theta_s_vector, phi_s_vector, surfParam,TRAp_Betac,sigmaParam_Betac);
% 

    pr(i) = pr_scalar;
    AbsE(i) = sqrt((pr(i) - P_IT)^2);
    beta_pr(i,:) = [beta_c(i), pr_scalar, AbsE(i)];
    coeff_matrix_all{i} = coeff_matrix;

end


% pr_db = db(pr,'power');
% P_IT_db = db(P_IT,'power');
% AbsE_db = sqrt((pr_db - P_IT_db).^2);
% loglog(beta_c,RMSE_db_2)

[~, idx_min] = min(AbsE);
beta_c_opt = beta_c(idx_min);

%% Save result - Temp

BetacResult.P_IT = P_IT;
BetacResult.PT = Pt;
BetacResult.beta_c_vs_pr = beta_pr; 
BetacResult.beta_c_opt = beta_c_opt;
BetacResult.coeff_matrix_all = coeff_matrix_all;
save('C:\Users\jalef\Downloads\beta2_RMSE\beam10\BetaC_H1000_Inc46.mat', 'BetacResult');
%%

%% plot
figure;
loglog(beta_c, AbsE, 'b-', 'LineWidth', 2);
xlabel('\beta_c (deg)');
ylabel('RMSE');
title('RMSE vs \beta_2 - \beta_0 = 0.03');
grid on;

%% plot

% 
%  radarEquation(R0, Rs, lambda, sigma_z, beta_c, Pt, D_rx, D_tx, ...
%                              Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
%                              SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
%                              beta_c, beta_c,theta_i, phi_i, theta_s_vector, phi_s_vector) ...
%                 - P_IT).^2, 'omitnan');

% opts = optimset('TolX',1e-4,'Display','off');

%% plot


% beta_c_range = 0.001:0.01:5;
% rmse_values = zeros(size(beta_c_range));
% 
% parfor i = 1:length(beta_c_range)
%     rmse_values(i) = rmseFun(beta_c_range(i));
% end
% 
% % Plot
% figure;
% plot(beta_c_range, rmse_values, 'b-', 'LineWidth', 2);
% xlabel('\beta_c (deg)');
% ylabel('RMSE');
% title('RMSE vs \beta_c');
% grid on;

%%

% beta_c_opt = fminbnd(rmseFun, 0.001, 20, opts);

end


function [Pr, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, beta_c, ...
                                   Pt, Ny_surf, Nx_surf, ...
                                   freq, SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
                                   beta0, beta2,theta_i, phi_i, theta_s_vector, phi_s_vector,surfParam,TRAp,sigmaParam_Betac)


% function Pr = radarEquation(R0, Rs, lambda, sigma_z, beta_c, ...
%                             theta_s, phi_s, Pt, D_rx, D_tx, Ny_surf, Nx_surf, ...
%                             freq, SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
%                             beta0, beta2, theta_i, phi_i, theta_s_vector , phi_s_vector)

  % Integrate coherent NRCS over the surface patch
  % Patch around specular point of size Ny_surf×Nx_surf
  x_vec  = surfParam.Xplot;
  y_vec  = surfParam.Yplot;
  [X,Y]  = meshgrid(x_vec, y_vec);

  % Ranges to each element
%   R1 = zeros(size(Nx_surf,Ny_surf))+sigmaParam_Betac.R0  ;
  R1   = repmat(sigmaParam_Betac.R0 ,Ny_surf,Nx_surf,1);

  R2 = sigmaParam_Betac.r2mat;
  R2_vector = R2(:);

  % Local angles
  theta1 = theta_i;
  theta2 = theta_s_vector;
  phi1   = 0;
  phi2   = phi_s_vector;

% theta_s_vector
% phi_s_vector

% theta_i=10; %deg
% phi_i=0;%deg
% theta_s=10; %deg
% phi_s=0;%deg
% f=1.41;  %GHz
% SM=15; %[5 20 40];
% sstd=1; %[0.5 1.5 3];  %cm
% R0=22000000./cosd(theta_i); %m
% R2=1500./cosd(theta_i); %m
% betac=1; %[deg] %half beamwidth
% beta0=betac;
% beta2=betac;
% coeff=coerente(f*10^9,SM,sstd/100,R0,R2,beta0rad,beta2rad,t,p,ts,ps);

% theta_s =10;
% phi_s = 2;
% % sigma_z=1;
% % theta_i = 10;
% beta_c = 10;
% beta0 = beta_c;
% beta2 = beta_c;
% theta_i = 10;
freq_hz = freq*10^9;
beta0_rad = deg2rad(beta0);
beta2_rad = deg2rad(beta2);
theta_i_rad = deg2rad(theta_i);
phi_i_rad = deg2rad(phi_i);

theta_s_vector_filtered = theta_s_vector;
theta_s_vector_filtered(theta_s_vector > 70) = NaN;

theta_s_rad = deg2rad(theta_s_vector_filtered);
phi_s_rad = deg2rad(phi_s_vector);

coeff_vector = NaN(size(theta_s_rad));

% g0 = HSAVERSTXgain(TRAp);
% gs = HSAVERSRXgain(TRAp,surfParam.Xplot,surfParam.Yplot,surfParam.Nx_surf,surfParam.Ny_surf);
% gs = gs(:);

for i = 1:length(R2_vector)
    coeff_vector(i) = coerente(freq_hz, SMvol, sigma_z, R0_slant, R2_vector(i), ...
                               beta0_rad, beta2_rad, theta_i_rad, phi_i_rad, ...
                               theta_s_rad(i), phi_s_rad(i));
end


% coeff_vector = abs(coeff_vector).^2;
coeff_matrix = reshape(coeff_vector, size(R2));

% coeff_matrix_power2 = abs(coeff_matrix).^2;
%test amir$
% coeff_matrix = abs(coeff_matrix).^2;
%%

coeff_matrix = abs(coeff_matrix);

% is_nan_complex = isnan(real(coeff_matrix)) & isnan(imag(coeff_matrix));
% coeff_matrix(is_nan_complex) = NaN;
% zero_imag = imag(coeff_matrix) == 0;
% coeff_matrix(zero_imag) = real(coeff_matrix(zero_imag));

% coeff_matrix = abs(coeff_matrix).^2;
%amirtest
% coeff_matrix_test = zeros(size(coeff_matrix))+db2pow(30.124);

theta_s_matrix = reshape(theta_s_vector, size(R2));
phi_s_matrix = reshape(phi_s_vector, size(R2));

  theta2 = theta_s_matrix;
  phi2   = phi_s_matrix;

 center_row = ceil(size(theta_s_matrix, 1) / 2);
center_col = ceil(size(theta_s_matrix, 2) / 2);
theta2_SP = theta_s_matrix(center_row, center_col);
theta1_SP = theta_i;

% filter high angles (change to Nan)
%check coerente works with vector
% coeff = coerente(freq_hz, SMvol, sigma_z_cm/100, R0_slant, Rs_slant, ...
%                  beta0_rad, beta2_rad, theta_i_rad, phi_i_rad, ...
%                  theta_s_rad, phi_s_rad,R0,Rs);

% sigma_c = abs(coeff)^2; %check with prof. (coeff for ph
% db(sigma_c,'power')

  dA        = (x_vec(2)-x_vec(1)) * (y_vec(2)-y_vec(1));
  integrand = HSAVERSTXgain(TRAp) .* HSAVERSRXgain(TRAp,surfParam.Xplot,surfParam.Yplot,surfParam.Nx_surf,surfParam.Ny_surf) ...
            .* coeff_matrix ./ (R1.^2 .* R2.^2); % .* sigma_c

  Pr = Pt * (lambda^2/(4*pi)^3) * sum(integrand(:),'omitnan') * abs(dA);
end

% function gain = getGainMatrix(GAIN_MATRIX, theta_mat, phi_mat)
%     theta_deg = max(0, min(90, theta_mat));
%     phi_deg   = mod(phi_mat, 360);
% 
%     inc_idx = round(theta_deg * 10) + 1;
%     azi_idx = round(phi_deg   * 10) + 1;
% 
%     linear_idx = sub2ind(size(GAIN_MATRIX), azi_idx, inc_idx);
%     gain = reshape(GAIN_MATRIX(linear_idx), size(theta_deg));
% 
% end

%% HSAVERS TX gain

function TxGain = HSAVERSTXgain(TRAp)

GPS_theta         = TRAp.theta_Txpattern;
GPS_phi           = TRAp.phi_Txpattern; 
GPSPattern_linear = TRAp.GainTx_linear(:,:,1);

TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles(3) +  TRAp.Yaw_TX;

[Txmat_ECEF2OF, Txmat_OF2BF, Txmat_BFtoAF] = matrixFrameGen(TRAp.TxECEF,  TRAp.VTxECEF, TRAp.Tx_nomAtt_angles, TRAp.TxAttitude_angles, TRAp.TxEuler_angles);

TRx_K = (TRAp.RxECEF - TRAp.TxECEF);

TRx_OF = Txmat_ECEF2OF*TRx_K;
TRx_BF = Txmat_OF2BF*TRx_OF;
TRx_AF = Txmat_BFtoAF*TRx_BF;

%piercing angles
thetaGPS_pier = acosd(dot(TRx_AF./norm(TRx_AF), [0,0,1]));
phiGPS_pier   = atan2d(TRx_AF(2), TRx_AF(1));

%DC: comment added on oct 2023
% [~,ind_phiTransm] = min(abs(abs(phiGPS_pier) - GPS_phi)); 
[~,ind_phiTransm] = min(abs(phiGPS_pier - GPS_phi));
[~,ind_thTransm]  = min(abs(thetaGPS_pier    - GPS_theta));

% UP_pattern     = antennaPattern_linear(ind_phiUP, ind_thUP);
TxGain = GPSPattern_linear(ind_phiTransm, ind_thTransm);

end


%% RX gain
function RxGain = HSAVERSRXgain(TRAp,Xplot,Yplot,Nx_surf,Ny_surf)
[regGridSim_x, regGridSim_y] = meshgrid(Xplot,Yplot);
regGridSim_z = zeros(Nx_surf,Ny_surf);
xp = regGridSim_x;
yp = regGridSim_y;
zp = regGridSim_z;

pSGR(:,:,1) = xp;
pSGR(:,:,2) = yp;
pSGR(:,:,3) = zp;     
        
pENU    = rotation_counterclockwise_mat(pSGR, TRAp.phi_ENU2local, 'rad');   
Rx_ENU  = rotation_counterclockwise(TRAp.Rxlocal, TRAp.phi_ENU2local, 'rad');   
VRx_ENU = rotation_counterclockwise(TRAp.VRxlocal, TRAp.phi_ENU2local, 'rad');  
            
[ph,  lm,  ~] = xyz2llh(TRAp.OENUinECEF(1), TRAp.OENUinECEF(2), TRAp.OENUinECEF(3));
            
Amat_ECEF2ENU = [ -sin(lm)                cos(lm)       0;
                  -sin(ph)*cos(lm)  -sin(ph)*sin(lm)  cos(ph);...
                   cos(ph)*cos(lm)   cos(ph)*sin(lm)  sin(ph)];
            
vdiff         =  [0,0,0].' -  TRAp.OENUinECEF;
OecefINenu    = Amat_ECEF2ENU*vdiff;
            
Amat_ENU2ECEF = Amat_ECEF2ENU.';
            
pdiffENU(:,:,1) = pENU(:,:,1) - OecefINenu(1);
pdiffENU(:,:,2) = pENU(:,:,2) - OecefINenu(2);
pdiffENU(:,:,3) = pENU(:,:,3) - OecefINenu(3);

pRxdiffENU      = Rx_ENU - OecefINenu;
         

pECEF(:,:,1) = -sin(lm).*pdiffENU(:,:,1) - sin(ph)*cos(lm).*pdiffENU(:,:,2) + cos(ph)*cos(lm).*pdiffENU(:,:,3);
pECEF(:,:,2) =  cos(lm).*pdiffENU(:,:,1) - sin(ph)*sin(lm).*pdiffENU(:,:,2) + cos(ph)*sin(lm).*pdiffENU(:,:,3);
pECEF(:,:,3) =            0              + cos(ph).*pdiffENU(:,:,2)         + sin(ph).*pdiffENU(:,:,3);  

RxECEF  = Amat_ENU2ECEF*pRxdiffENU;
vRxECEF = Amat_ENU2ECEF*VRx_ENU;
[Amat_ECEF2OF, Amat_OF2BF, Amat_BFtoAF] = matrixFrameGen(RxECEF, vRxECEF, TRAp.Rx_nomAtt_angles, TRAp.RxAttitude_angles, TRAp.RxEuler_angles);
        
rs_ECEF(:,:,1) = (pECEF(:,:,1) - RxECEF(1));
rs_ECEF(:,:,2) = (pECEF(:,:,2) - RxECEF(2));
rs_ECEF(:,:,3) = (pECEF(:,:,3) - RxECEF(3));
        
rs_OF(:,:,1) = Amat_ECEF2OF(1,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(1,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(1,3).*rs_ECEF(:,:,3);
rs_OF(:,:,2) = Amat_ECEF2OF(2,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(2,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(2,3).*rs_ECEF(:,:,3);
rs_OF(:,:,3) = Amat_ECEF2OF(3,1).*rs_ECEF(:,:,1) + Amat_ECEF2OF(3,2).*rs_ECEF(:,:,2) + Amat_ECEF2OF(3,3).*rs_ECEF(:,:,3);

rs_BF(:,:,1) = Amat_OF2BF(1,1).*rs_OF(:,:,1) + Amat_OF2BF(1,2).*rs_OF(:,:,2) + Amat_OF2BF(1,3).*rs_OF(:,:,3);
rs_BF(:,:,2) = Amat_OF2BF(2,1).*rs_OF(:,:,1) + Amat_OF2BF(2,2).*rs_OF(:,:,2) + Amat_OF2BF(2,3).*rs_OF(:,:,3);
rs_BF(:,:,3) = Amat_OF2BF(3,1).*rs_OF(:,:,1) + Amat_OF2BF(3,2).*rs_OF(:,:,2) + Amat_OF2BF(3,3).*rs_OF(:,:,3);

rs_AF(:,:,1) = Amat_BFtoAF(1,1).*rs_BF(:,:,1) + Amat_BFtoAF(1,2).*rs_BF(:,:,2) + Amat_BFtoAF(1,3).*rs_BF(:,:,3);
rs_AF(:,:,2) = Amat_BFtoAF(2,1).*rs_BF(:,:,1) + Amat_BFtoAF(2,2).*rs_BF(:,:,2) + Amat_BFtoAF(2,3).*rs_BF(:,:,3);
rs_AF(:,:,3) = Amat_BFtoAF(3,1).*rs_BF(:,:,1) + Amat_BFtoAF(3,2).*rs_BF(:,:,2) + Amat_BFtoAF(3,3).*rs_BF(:,:,3);
                   
thetas_el_pier   = asind(rs_AF(:,:,3)./sqrt(dot(rs_AF,rs_AF,3)));

thetas_zen_pier  = acosd(rs_AF(:,:,3)./sqrt(dot(rs_AF,rs_AF,3))); 

phis_pier       = atan2d(rs_AF(:,:,2), rs_AF(:,:,1));

phiRxDown_3D      = reshape(TRAp.phiRxDown,1,1,length(TRAp.phiRxDown));
phiRxDown_array   = repmat(phiRxDown_3D,Ny_surf,Nx_surf,1);
          
thetaRxDown_3D    = reshape(TRAp.thetaRxDown,1,1,length(TRAp.thetaRxDown));
thetaRxDown_array = repmat(thetaRxDown_3D,Ny_surf,Nx_surf,1);

[~, phi_ind] = min(abs(phiRxDown_array    - phis_pier),[],3);
[~, the_ind] = min(abs(thetaRxDown_array  - thetas_zen_pier),[],3);        
indG1D       = sub2ind([size(TRAp.GainRxDown_LHCP,1), size(TRAp.GainRxDown_LHCP,2)], reshape(phi_ind, Nx_surf*Ny_surf, 1), reshape(the_ind, Nx_surf*Ny_surf, 1)); 

%%
linear_idx = sub2ind(size(TRAp.GainRxDown_LHCP), phi_ind, the_ind);
RxGain = reshape(TRAp.GainRxDown_LHCP(linear_idx), size(phi_ind));

end
