function beta_opt = findBeta0Beta2_New(R0, Rs, lambda, Pt, sigma_z, ...
                                       TRApGain_rx, Ny_surf, Nx_surf, TRApGain_tx, ...
                                       sigmaParam_Betac, TRAp_Betac, SMvol, ...
                                       gammaParam, surfParam)

Pt = db2pow(Pt);  % dB to linear
sigma_z_cm = sigma_z;
sigma_z = sigma_z / 100;

theta_i = TRAp_Betac.incAngle_enuRot;
phi_i = 0;

theta_s = sigmaParam_Betac.TH_s;
phi_s = sigmaParam_Betac.PH_s;

R0_slant = R0 ./ cosd(theta_i);
Rs_slant = Rs ./ cosd(theta_i);

%% Fresnel coefficient
cd = costantedielettrica(SMvol./100, TRAp_Betac.TxFrequency);
refr_index = sqrt(cd);
[rH, rV] = fresnel_Matlab(1, refr_index, theta_i);
reflectioncoeff = (rV + rH) / 2;

P_IT = Pt * gammaParam.Gain_LHCP_onsurf_lin_atSP * TRAp_Betac.TxGainAtSP ...
       * lambda^2 * abs(reflectioncoeff)^2 ...
       / ((4*pi)^2 * (R0_slant + Rs_slant)^2);

%% Flatten theta_s and phi_s for vectorised access
theta_s_vector = theta_s(:);
phi_s_vector = phi_s(:);
N = numel(theta_s_vector);

%% Grid of beta0 and beta2
beta0_vals = 0.001:2:100;
beta2_vals = 0.001:2:100;

RMSE_matrix = zeros(length(beta0_vals), length(beta2_vals));
Pr_matrix = zeros(length(beta0_vals), length(beta2_vals));

coeff_matrix_all = cell(length(beta0_vals), length(beta2_vals));

for i = 1:length(beta0_vals)
    for j = 1:length(beta2_vals)
        [pr_scalar, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, ...
            [], Pt, Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
            SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
            beta0_vals(i), beta2_vals(j), theta_i, phi_i, ...
            theta_s_vector, phi_s_vector, surfParam, TRAp_Betac, sigmaParam_Betac);

        Pr_matrix(i, j) = pr_scalar;
        RMSE_matrix(i, j) = sqrt((pr_scalar - P_IT)^2);
        coeff_matrix_all{i, j} = coeff_matrix;
    end
end

%% Find minimum RMSE
[min_val, linear_idx] = min(RMSE_matrix(:));
[idx_i, idx_j] = ind2sub(size(RMSE_matrix), linear_idx);
beta0_opt = beta0_vals(idx_i);
beta2_opt = beta2_vals(idx_j);
% Pr_matrix(idx_i,idx_j);
beta_opt.beta0 = beta0_opt;
beta_opt.beta2 = beta2_opt;
beta_opt.P_IT = P_IT;
beta_opt.RMSE_matrix = RMSE_matrix;
beta_opt.beta0_vals = beta0_vals;
beta_opt.beta2_vals = beta2_vals;
beta_opt.coeff_matrix_all = coeff_matrix_all;

% %% Save result
% save('C:\Users\jalef\Downloads\beta_C_RMSE\varnew3\Beta0_Beta2_Result.mat', 'beta_opt');
% 
% %% Plot
% figure;
% imagesc(beta2_vals, beta0_vals, RMSE_matrix);
% set(gca, 'YDir', 'normal');
% xlabel('\beta_2 (deg)');
% ylabel('\beta_0 (deg)');
% title('RMSE between P_r and P_{IT}');
% colorbar;
end


%% RADAR Equation

function [Pr, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, beta_c, ...
                                   Pt, Ny_surf, Nx_surf, ...
                                   freq, SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
                                   beta0, beta2,theta_i, phi_i, theta_s_vector, phi_s_vector,surfParam,TRAp,sigmaParam_Betac)

  x_vec  = surfParam.Xplot;
  y_vec  = surfParam.Yplot;
  [X,Y]  = meshgrid(x_vec, y_vec);

  R1   = repmat(sigmaParam_Betac.R0 ,Ny_surf,Nx_surf,1);

  R2 = sigmaParam_Betac.r2mat;
  R2_vector = R2(:);

  % Local angles
  theta1 = theta_i;
  theta2 = theta_s_vector;
  phi1   = 0;
  phi2   = phi_s_vector;

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

coeff_matrix = reshape(coeff_vector, size(R2));

coeff_matrix = abs(coeff_matrix);

theta_s_matrix = reshape(theta_s_vector, size(R2));
phi_s_matrix = reshape(phi_s_vector, size(R2));

  theta2 = theta_s_matrix;
  phi2   = phi_s_matrix;

 center_row = ceil(size(theta_s_matrix, 1) / 2);
center_col = ceil(size(theta_s_matrix, 2) / 2);
theta2_SP = theta_s_matrix(center_row, center_col);
theta1_SP = theta_i;

  dA        = (x_vec(2)-x_vec(1)) * (y_vec(2)-y_vec(1));
  integrand = HSAVERSTXgain(TRAp) .* HSAVERSRXgain(TRAp,surfParam.Xplot,surfParam.Yplot,surfParam.Nx_surf,surfParam.Ny_surf) ...
            .* coeff_matrix ./ (R1.^2 .* R2.^2); % .* sigma_c

  Pr = Pt * (lambda^2/(4*pi)^3) * sum(integrand(:),'omitnan') * abs(dA);
end


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
