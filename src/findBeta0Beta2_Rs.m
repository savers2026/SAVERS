function beta_optx = findBeta0Beta2_Rs(R0, Rs, lambda, Pt, sigma_z, ...
                                       TRApGain_rx, Ny_surf, Nx_surf, TRApGain_tx, ...
                                       sigmaParam_Betac, TRAp_Betac, SMvol, ...
                                       gammaParam, surfParam,geoSYSp)

Pt = db2pow(Pt);  % dB to linear
sigma_z = 0;
sigma_z_cm = sigma_z;
sigma_z = sigma_z / 100;



beta0_val = 0.03;
beta2_vals = 0.001:0.1:100;
theta_i_org = TRAp_Betac.incAngle_enuRot;
theta_i_range = 5:10:65;
phi_i = 0;


%% Rs org in center
Rs_org = Rs;
% num_each_side = 2;  % values below and above Rs
% num_total_samples = 2 * num_each_side + 1;
% 
% range_min = Rs_org / 3;
% range_max = Rs_org * 3;
% 
% step_down = (Rs_org - range_min) / num_each_side;
% step_up   = (range_max - Rs_org) / num_each_side;
% 
% lower_values = Rs_org - step_down * (num_each_side:-1:1);
% upper_values = Rs_org + step_up   * (1:num_each_side);
% 
% Rs_base_range = [lower_values, Rs_org, upper_values];
% 
% Rs_base_range_rounded = round(Rs_base_range, -2);
% Rs_range = unique(Rs_base_range_rounded);


% Rs_range = [500 , 1000, 2000 , 3000];
% Rs_range = Rs_org; % to disable iteration !

num_beta0  = numel(beta0_val);
num_beta2  = numel(beta2_vals);
num_theta     = numel(theta_i_range);
ABSe_cube  = zeros(num_theta , num_beta2 , num_beta0);
Pr_Theta_beta2_matrix = zeros(num_theta, num_beta2);

for Theta_idx = 1:length(theta_i_range)

theta_i = theta_i_range(Theta_idx);
Rs = Rs_org;

R0_slant = R0 ./ cosd(theta_i);
Rs_slant = Rs ./ cosd(theta_i);


Xmax_Range =  round(sqrt(2*2.99792458e2*(geoSYSp.Dly_RNGmax)*Rs*cosd(theta_i))/cosd(theta_i)^2/1000)*1000;
Xmin_Range   = -Xmax_Range; 
XplotPos = 0.0:geoSYSp.dx:Xmax_Range;
XplotNeg = -(0.0+geoSYSp.dx:geoSYSp.dx:abs(Xmin_Range));
Xplot    = [fliplr(XplotNeg) XplotPos];
Nx_surf  = length(Xplot); 

Yabs_Range   =  Xmax_Range;
YplotPos = 0.0:geoSYSp.dy:Yabs_Range;
YplotNeg = -YplotPos(2:length(YplotPos));
Yplot    = [fliplr(YplotNeg) YplotPos];
Yplot    = flipud(Yplot');
Ny_surf  = length(Yplot);


%% =====

  x_vec  = Xplot;
  y_vec  = Yplot;
  [X,Y]  = meshgrid(x_vec, y_vec);

[~,~,R1_scalar] = xyz2llh(TRAp_Betac.TxECEF(1),TRAp_Betac.TxECEF(2),TRAp_Betac.TxECEF(3));

  R1   = repmat(R1_scalar ,Ny_surf,Nx_surf,1);

%% TEST Height

[regGridSim_x, regGridSim_y] = meshgrid(x_vec,y_vec);
regGridSim_z = zeros(Nx_surf,Ny_surf);
xp = regGridSim_x;
yp = regGridSim_y;
zp = regGridSim_z;
clear pSGR RxPdiff
pSGR(:,:,1) = xp;
pSGR(:,:,2) = yp;
pSGR(:,:,3) = zp;      



%%


[lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, Height_Rx_Temp_Rs] = xyz2llh(TRAp_Betac.RxECEF(1),TRAp_Betac.RxECEF(2),TRAp_Betac.RxECEF(3));
lat_Rx_Temp_Rs = rad2deg(lat_Rx_Temp_Rs);
lon_Rx_Temp_Rs = rad2deg(lon_Rx_Temp_Rs);
RxLLA_Rs_temp = [lat_Rx_Temp_Rs,lon_Rx_Temp_Rs,Rs];
[RxECEF_Rs] = lla2ecef(RxLLA_Rs_temp)';
%
[OENUinLLA] = Specular_Point_Position_Calculation(RxECEF_Rs', TRAp_Betac.TxECEF');
OENUinECEF = lla2ecef(OENUinLLA);
OENUinECEF = OENUinECEF';
%
TxENU_Rs  = ecef2enu_dc(TRAp_Betac.TxECEF, OENUinECEF); 
RxENU_Rs  = ecef2enu_dc(RxECEF_Rs, OENUinECEF);
phiTx = atan2(TxENU_Rs(2), TxENU_Rs(1));
phi_ENU2local = pi - phiTx;
Txlocal_rs  = rotation_clockwise(TxENU_Rs,   phi_ENU2local, 'rad');
Rxlocal_rs  = rotation_clockwise(RxENU_Rs,   phi_ENU2local, 'rad');

Cords.Txlocal_rs = Txlocal_rs;
Cords.Rxlocal_rs = Rxlocal_rs;
Cords.OENUinECEF = OENUinECEF;
Cords.phi_ENU2local = phi_ENU2local;
Cords.TxENU_Rs = TxENU_Rs;
Cords.RxENU_Rs = RxENU_Rs;
Cords.RxECEF_Rs = RxECEF_Rs;


%%

RxPdiff(:,:,1) = Cords.Rxlocal_rs(1) - pSGR(:,:,1); 
RxPdiff(:,:,2) = Cords.Rxlocal_rs(2) - pSGR(:,:,2); 
RxPdiff(:,:,3) = Cords.Rxlocal_rs(3) - pSGR(:,:,3); 
RxPmod  = sqrt(RxPdiff(:,:,1).^2 + RxPdiff(:,:,2).^2 + RxPdiff(:,:,3).^2);    

theta_s = asind(sqrt((pSGR(:,:,1) - Cords.Rxlocal_rs(1)).^2 + (pSGR(:,:,2) - Cords.Rxlocal_rs(2)).^2)./RxPmod);
phi_s = atan2d(Cords.Rxlocal_rs(2)  - pSGR(:,:,2), Cords.Rxlocal_rs(1) - pSGR(:,:,1));  

[RxGain,Gain_LHCP_onsurf_lin_atSP,Gain_RHCP_onsurf_lin_atSP] = HSAVERSRXgain(TRAp_Betac,Xplot,Yplot,Nx_surf,Ny_surf,Rs,Cords);

Rx_Gain_All.RxGain                     = RxGain;
Rx_Gain_All.Gain_LHCP_onsurf_lin_atSP  = Gain_LHCP_onsurf_lin_atSP;
Rx_Gain_All.Gain_RHCP_onsurf_lin_atSP  = Gain_RHCP_onsurf_lin_atSP;

TxGain_Rs = HSAVERSTXgain(TRAp_Betac,Rs,Cords);

%% Fresnel coefficient
cd = costantedielettrica(SMvol./100, TRAp_Betac.TxFrequency);
refr_index = sqrt(cd);
[rH, rV] = fresnel_Matlab(1, refr_index, theta_i);
reflectioncoeff = (rV + rH) / 2;


% Note: Tx Gain on SP doest not relevant to Rx height, So no change on TxGainAtSP

P_IT = Pt * Gain_LHCP_onsurf_lin_atSP * TxGain_Rs ...
       * lambda^2 * abs(reflectioncoeff)^2 ...
       / ((4*pi)^2 * (R0_slant + Rs_slant)^2);


%% Flatten theta_s and phi_s for vectorised access
theta_s_vector = theta_s(:);
phi_s_vector = phi_s(:);
% N = numel(theta_s_vector);

%% Beta2 Iteration - Fixing Beta0 0.03 deg



ABSe_matrix = zeros(length(beta2_vals), length(beta0_val));
Pr_matrix = zeros(length(beta2_vals), length(beta0_val));
coeff_matrix_all = cell(length(beta2_vals), length(beta0_val));

for i = 1:length(beta0_val)
    for j = 1:length(beta2_vals)
        [pr_scalar, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, ...
            [], Pt, Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
            SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
            beta0_val(i), beta2_vals(j), theta_i, phi_i, ...
            theta_s_vector, phi_s_vector, surfParam, TRAp_Betac, sigmaParam_Betac,geoSYSp,Rx_Gain_All,RxPmod,Cords);

        Pr_matrix(j, i) = pr_scalar;
        ABSe_matrix(j, i) = abs(pr_scalar - P_IT);
        coeff_matrix_all{j, i} = coeff_matrix;
        beta_pr(j,:) = [beta2_vals(j), pr_scalar, ABSe_matrix(j)];
    end
end

P_IT_matrix = repmat(P_IT, size(Pr_matrix));

% Compute RMSE per Beta0 (columns) and per Beta2 (rows)
RMSE_matrix_col_beta0 = sqrt(mean((Pr_matrix - P_IT_matrix).^2, 1));  % Beta0 along columns
RMSE_matrix_row_beta2 = sqrt(mean((Pr_matrix - P_IT_matrix).^2, 2));  % Beta2 along rows

%  tolerance based selection for Beta0 - Avoiding floating point error
min_col_val = min(RMSE_matrix_col_beta0);
tolerance_col = 1e-4 * min_col_val;
close_col_indices = find(abs(RMSE_matrix_col_beta0 - min_col_val) <= tolerance_col);
if isempty(close_col_indices)
    [~, idxbeta0] = min(RMSE_matrix_col_beta0);
else
    idxbeta0 = close_col_indices(1);
end

%  tolerance based selection for Beta2 - Avoiding floating point error
min_row_val = min(RMSE_matrix_row_beta2);
tolerance_row = 1e-4 * min_row_val;
close_row_indices = find(abs(RMSE_matrix_row_beta2 - min_row_val) <= tolerance_row);
if isempty(close_row_indices)
    [~, idxbeta2] = min(RMSE_matrix_row_beta2);
else
    idxbeta2 = close_row_indices(1);
end

optimal_beta0 = beta0_val(idxbeta0);
optimal_beta2 = beta2_vals(idxbeta2);

% plot5

% figure;
% imagesc(beta0_val, beta2_vals, db(ABSe_matrix,'power'));
% % caxis([4.9956e-21, 1.4185e-14]); %temp
% xlabel('\beta_0');
% ylabel('\beta_2');
% title(sprintf('ABSe vs. \\beta_0 and \\beta_2  |  Height (m): %.2f', Rs));
% colorbar;
% axis xy; 
% 
% hold on;
% plot(optimal_beta0, optimal_beta2, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
% text(optimal_beta0, optimal_beta2, sprintf('  (%.2f, %.2f)', optimal_beta0, optimal_beta2), ...
%     'Color', 'black', 'FontSize', 10, 'FontWeight', 'bold');

Pr_Theta_beta2_matrix(Theta_idx, :) = Pr_matrix(:, 1).';
optimal_beta2_all(Theta_idx, :) = optimal_beta2;
P_IT_matrix_all(Theta_idx, :) = P_IT_matrix';
ABSe_cube(Theta_idx,:,:) = ABSe_matrix;
end

%% END Beta2 iteration


RMSE_beta2 = squeeze( sqrt( mean( ABSe_cube.^2 , 1 ) ) );
[~,idx_best] = min(RMSE_beta2);
best_beta2 = beta2_vals(idx_best);

% %% AVG beta2 is the best?
% best_beta2 = mean(optimal_beta2_all);
% %%
% mean_row_Pr_Theta_beta2_matrix = mean(Pr_Theta_beta2_matrix,1);
% 
% % discard very different beta2, logic is if the beta2 is greater than twice of mean, then discard
% optimal_beta2_mean = mean(optimal_beta2_all);
% % beta2_discard_idx = find(optimal_beta2_all > 2 * optimal_beta2_mean);
% beta2_discard_idx = find(optimal_beta2_all > 100000000);
% 
% Pr_beta2_filtered = Pr_Theta_beta2_matrix;
% Pr_beta2_filtered(beta2_discard_idx, :) = [];
% P_IT_matrix_filtered = P_IT_matrix_all;
% P_IT_matrix_filtered(beta2_discard_idx, :) = [];
% 
% 
% mean_Pr_beta2 = mean(Pr_beta2_filtered, 1);  % mean over theta_i
% mean_P_IT = mean(P_IT_matrix_filtered, 1);  % mean over theta_i
% 
% RMSE_beta2 = sqrt((mean_Pr_beta2 - mean_P_IT).^2);  % RMSE per beta2
% [min_RMSE, idx_min_beta2] = min(RMSE_beta2);
% best_beta2 = beta2_vals(idx_min_beta2);
% 
% 
% RMSE_beta2_db = sqrt((db(mean_Pr_beta2,'power') - db(mean_P_IT,'power')).^2);  % RMSE per beta2
% [min_RMSE_db, idx_min_beta2_db] = min(RMSE_beta2_db);
% best_beta2_db = beta2_vals(idx_min_beta2_db);


%% ===========



optimal_beta0 = 0.03;
beta0_opt = optimal_beta0;
beta2_opt = best_beta2;
beta_c_opt = (beta0_opt+beta2_opt)/2;
beta_optx = [beta0_opt, beta2_opt, beta_c_opt];

%% final plot

pr_matrix_optBeta = db(Pr_Theta_beta2_matrix(:,idx_best), 'power');
P_IT_optBeta = db(P_IT_matrix_all(:,idx_best), 'power');

figure('Color', 'w'); hold on; grid on;

plot(theta_i_range, P_IT_optBeta, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'P_{IT}');
plot(theta_i_range, pr_matrix_optBeta, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', sprintf('p_r (at \\beta_2 = %.2f°)', beta2_vals(idx_best)));

xlabel('Height (m)', 'FontSize', 10);
ylabel('Power (dB)', 'FontSize', 14);
ylim([min(min(pr_matrix_optBeta),min(P_IT_optBeta))-1, max(max(pr_matrix_optBeta),max(P_IT_optBeta))+1]);
title(sprintf('p_r and P_{IT} vs. Incidence angle at Optimal \\beta_2 (H: %d m)', round(Rs_org)), ...
    'FontSize', 10);
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10);

%% final plot


% % Save results
% beta_opt.beta_c_vs_pr = beta_pr;
% beta_opt.beta0 = beta0_opt;
% beta_opt.beta2 = optimal_beta2;
% beta_opt.betac = beta_c_opt;
% beta_opt.P_IT = P_IT;
% beta_opt.RMSE_matrix = RMSE_matrix_row_beta2;
% beta_opt.beta0_vals = beta0_val;
% beta_opt.beta2_vals = beta2_vals;
% beta_opt.coeff_matrix_all = coeff_matrix_all;


%% Save result
% save('C:\Users\jalef\Downloads\betc_new_july\h4\BetaC_H500_Inc16.mat', 'beta_opt');

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
                                   beta0, beta2,theta_i, phi_i, theta_s_vector, phi_s_vector,surfParam,TRAp,sigmaParam_Betac,geoSYSp,Rx_Gain_All,RxPmod,Cords)

%% Dynamic area

Xmax_Range =  round(sqrt(2*2.99792458e2*(geoSYSp.Dly_RNGmax)*Rs*cosd(theta_i))/cosd(theta_i)^2/1000)*1000;
Xmin_Range   = -Xmax_Range; 
XplotPos = 0.0:geoSYSp.dx:Xmax_Range;
XplotNeg = -(0.0+geoSYSp.dx:geoSYSp.dx:abs(Xmin_Range));
Xplot    = [fliplr(XplotNeg) XplotPos];
Nx_surf  = length(Xplot); 

Yabs_Range   =  Xmax_Range;
YplotPos = 0.0:geoSYSp.dy:Yabs_Range;
YplotNeg = -YplotPos(2:length(YplotPos));
Yplot    = [fliplr(YplotNeg) YplotPos];
Yplot    = flipud(Yplot');
Ny_surf  = length(Yplot);


%% =====

  x_vec  = Xplot;
  y_vec  = Yplot;
  [X,Y]  = meshgrid(x_vec, y_vec);

  [~,~,R1_scalar] = xyz2llh(TRAp.TxECEF(1),TRAp.TxECEF(2),TRAp.TxECEF(3));
  R1   = repmat(R1_scalar ,Ny_surf,Nx_surf,1);

% %% TEST Height
% 
% [regGridSim_x, regGridSim_y] = meshgrid(x_vec,y_vec);
% regGridSim_z = zeros(Nx_surf,Ny_surf);
% xp = regGridSim_x;
% yp = regGridSim_y;
% zp = regGridSim_z;
% pSGR(:,:,1) = xp;
% pSGR(:,:,2) = yp;
% pSGR(:,:,3) = zp;      
% RxPdiff(:,:,1) = TRAp.Rxlocal(1) - pSGR(:,:,1); 
% RxPdiff(:,:,2) = TRAp.Rxlocal(2) - pSGR(:,:,2); 
% RxPdiff(:,:,3) = Rs - pSGR(:,:,3); 
% RxPmod  = sqrt(RxPdiff(:,:,1).^2 + RxPdiff(:,:,2).^2 + RxPdiff(:,:,3).^2);    

%% TEST Height


%   R2 = sigmaParam_Betac.r2mat;
  R2 = RxPmod;
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
  integrand = HSAVERSTXgain(TRAp,Rs,Cords) .* Rx_Gain_All.RxGain ...
            .* coeff_matrix ./ (R1.^2 .* R2.^2); % .* sigma_c

  Pr = ((Pt * lambda^2)/(4*pi)^3) * sum(integrand(:),'omitnan') * abs(dA);
end


%% HSAVERS TX gain

function TxGain = HSAVERSTXgain(TRAp,Rs,Cords)

GPS_theta         = TRAp.theta_Txpattern;
GPS_phi           = TRAp.phi_Txpattern; 
GPSPattern_linear = TRAp.GainTx_linear(:,:,1);

TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles(3) +  TRAp.Yaw_TX;

[Txmat_ECEF2OF, Txmat_OF2BF, Txmat_BFtoAF] = matrixFrameGen(TRAp.TxECEF,  TRAp.VTxECEF, TRAp.Tx_nomAtt_angles, TRAp.TxAttitude_angles, TRAp.TxEuler_angles);

% Added for updating Rs
% [lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, Height_Rx_Temp_Rs] = xyz2llh(TRAp.RxECEF(1),TRAp.RxECEF(2),TRAp.RxECEF(3));
% lat_Rx_Temp_Rs = rad2deg(lat_Rx_Temp_Rs);
% lon_Rx_Temp_Rs = rad2deg(lon_Rx_Temp_Rs);
% RxLLA_Rs_temp = [lat_Rx_Temp_Rs,lon_Rx_Temp_Rs,Rs];
% [RxECEF_Rs] = lla2ecef(RxLLA_Rs_temp)';
%
RxECEF_Rs = Cords.RxECEF_Rs;

TRx_K = (RxECEF_Rs - TRAp.TxECEF);

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
function [RxGain,Gain_LHCP_onsurf_lin_atSP,Gain_RHCP_onsurf_lin_atSP] = HSAVERSRXgain(TRAp,Xplot,Yplot,Nx_surf,Ny_surf,Rs,Cords)
[regGridSim_x, regGridSim_y] = meshgrid(Xplot,Yplot);
regGridSim_z = zeros(Nx_surf,Ny_surf);
xp = regGridSim_x;
yp = regGridSim_y;
zp = regGridSim_z;

pSGR(:,:,1) = xp;
pSGR(:,:,2) = yp;
pSGR(:,:,3) = zp;     
OENUinECEF = Cords.OENUinECEF;

pENU    = rotation_counterclockwise_mat(pSGR, Cords.phi_ENU2local, 'rad');   
Rx_ENU  = rotation_counterclockwise(Cords.Rxlocal_rs, Cords.phi_ENU2local, 'rad');   
VRx_ENU = rotation_counterclockwise(TRAp.VRxlocal, Cords.phi_ENU2local, 'rad');  
            
[ph,  lm,  ~] = xyz2llh(OENUinECEF(1), OENUinECEF(2), OENUinECEF(3));
            
Amat_ECEF2ENU = [ -sin(lm)                cos(lm)       0;
                  -sin(ph)*cos(lm)  -sin(ph)*sin(lm)  cos(ph);...
                   cos(ph)*cos(lm)   cos(ph)*sin(lm)  sin(ph)];
            
vdiff         =  [0,0,0].' -  OENUinECEF;
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

% Gain on SP
Gain_LHCP_1D         = reshape(TRAp.GainRxDown_LHCP, size(TRAp.GainRxDown_LHCP,1).*size(TRAp.GainRxDown_LHCP,2),1);
Gain_LHCP_onsurf1D   = Gain_LHCP_1D(indG1D);
Gain_LHCP_onsurf_lin = reshape(Gain_LHCP_onsurf1D, Ny_surf, Nx_surf);

Gain_RHCP_1D         = reshape(TRAp.GainRxDown_RHCP, size(TRAp.GainRxDown_RHCP,1).*size(TRAp.GainRxDown_RHCP,2),1);
Gain_RHCP_onsurf1D   = Gain_RHCP_1D(indG1D);
Gain_RHCP_onsurf_lin = reshape(Gain_RHCP_onsurf1D, Ny_surf, Nx_surf);

indx_origin = Xplot == 0.0;
indy_origin = Yplot == 0.0;

Gain_LHCP_onsurf_lin_atSP = Gain_LHCP_onsurf_lin(indy_origin, indx_origin);
Gain_RHCP_onsurf_lin_atSP = Gain_RHCP_onsurf_lin(indy_origin, indx_origin);


end
