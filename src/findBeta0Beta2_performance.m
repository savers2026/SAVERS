function beta_optx = findBeta0Beta2_performance(R0, Rs, lambda, Pt, sigma_z, ...
    TRApGain_rx, Ny_surf, Nx_surf, TRApGain_tx, ...
    sigmaParam_Betac, TRAp_Betac, SMvol, ...
    gammaParam, surfParam, geoSYSp)

Pt = db2pow(Pt);
sigma_z_cm = 0;
sigma_z = 0;

beta0_val = 0.03;
beta2_vals = 0.1:0.1:80;
theta_i_range = 5:10:65;
% theta_i_range = 65;
phi_i = 0;
Rs_range = Rs;
Rs_org = Rs;
% Rs_range = [500,1000,2000,3000];

num_beta2 = numel(beta2_vals);
num_theta = numel(theta_i_range);
num_rs = numel(Rs_range);
ABSe_cube = zeros(num_theta, num_beta2, 1);
Pr_Theta_beta2_matrix = zeros(num_theta, num_beta2);



for Rs_idx = 1:num_rs
for Theta_idx = 1:num_theta
    theta_i = theta_i_range(Theta_idx);
    Rs = Rs_range(Rs_idx);
    R0_slant = R0 / cosd(theta_i);
    Rs_slant = Rs / cosd(theta_i);

    [x_vec, y_vec, X, Y, Nx_surf, Ny_surf, Area] = precomputeGrid(geoSYSp, Rs, TRAp_Betac,theta_i);

    clear pSGR
    [Cords, pSGR] = computeCords(Rs, TRAp_Betac, geoSYSp, Area,Rs_range,theta_i);
    clear RxPdiff
    RxPdiff(:,:,1) = Cords.Rxlocal_rs(1) - pSGR(:,:,1); 
    RxPdiff(:,:,2) = Cords.Rxlocal_rs(2) - pSGR(:,:,2); 
    RxPdiff(:,:,3) = Cords.Rxlocal_rs(3) - pSGR(:,:,3); 
    RxPmod  = sqrt(RxPdiff(:,:,1).^2 + RxPdiff(:,:,2).^2 + RxPdiff(:,:,3).^2);

    theta_s = asind(sqrt((pSGR(:,:,1) - Cords.Rxlocal_rs(1)).^2 + ...
                         (pSGR(:,:,2) - Cords.Rxlocal_rs(2)).^2) ./ RxPmod);
    phi_s = atan2d(Cords.Rxlocal_rs(2) - pSGR(:,:,2), ...
                   Cords.Rxlocal_rs(1) - pSGR(:,:,1));
    theta_s_vector = theta_s(:);
    phi_s_vector = phi_s(:);

    [RxGain, LHCP, RHCP] = HSAVERSRXgain(TRAp_Betac, Area.Xplot, Area.Yplot, Nx_surf, Ny_surf, Rs, Cords, Area,geoSYSp);
    Rx_Gain_All = struct('RxGain', RxGain, ...
                         'Gain_LHCP_onsurf_lin_atSP', LHCP, ...
                         'Gain_RHCP_onsurf_lin_atSP', RHCP);

    TxGain_Rs = HSAVERSTXgain(TRAp_Betac, Rs, Cords);
    cd = costantedielettrica(SMvol / 100, TRAp_Betac.TxFrequency);
    refr_index = sqrt(cd);
    [rH, rV] = fresnel_Matlab(1, refr_index, theta_i);
    reflectioncoeff = (rV + rH) / 2;

    P_IT = Pt * LHCP * TxGain_Rs * lambda^2 * abs(reflectioncoeff)^2 ...
           / ((4 * pi)^2 * (R0_slant + Rs_slant)^2 );
    
    Pr_matrix = zeros(num_beta2, 1);
    ABSe_matrix = zeros(num_beta2, 1);
    P_IT_matrix = repmat(P_IT, size(Pr_matrix));

    parfor j = 1:num_beta2
        [pr_scalar, ~] = radarEquation(R0, Rs, lambda, sigma_z, [], Pt, ...
            Ny_surf, Nx_surf, TRAp_Betac.TxFrequency, ...
            SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
            beta0_val, beta2_vals(j), theta_i, phi_i, ...
            theta_s_vector, phi_s_vector, ...
            surfParam, TRAp_Betac, sigmaParam_Betac, geoSYSp, ...
            Rx_Gain_All, RxPmod, Cords, Area);
        Pr_matrix(j) = pr_scalar;
        ABSe_matrix(j) = abs(pr_scalar - P_IT);
    end

    Pr_Theta_beta2_matrix(Theta_idx, :) = Pr_matrix;
    ABSe_cube(Theta_idx, :, 1) = ABSe_matrix;
    P_IT_matrix_all(Theta_idx, :) = P_IT_matrix';

%     figure;
%     plot(beta2_vals, db(Pr_Theta_beta2_matrix(Theta_idx,:),'power'), '-o', 'LineWidth', 2, 'MarkerSize', 3);
%     hold on;
%     plot(beta2_vals, db(P_IT_matrix_all(Theta_idx,:),'power'), '-s', 'LineWidth', 2, 'MarkerSize', 3);
%     legend({'Pr', 'P\_IT'}, 'Location', 'best');
%     xlabel('Beta2 (deg)');
%     ylabel('Power (dB)');
%     title(['Incidence angle (' num2str(theta_i_range(Theta_idx)) '°)']);
%     grid on;

end

RMSE_beta2 = squeeze(sqrt(mean(ABSe_cube.^2, 1)));
[~, idx_best_RMSE] = min(RMSE_beta2);
beta_optx = [beta0_val, beta2_vals(idx_best_RMSE)];

% Final plot
figure('Color', 'w'); hold on; grid on;
plot(theta_i_range, db(P_IT_matrix_all(:, idx_best_RMSE), 'power'), 'b-o', 'LineWidth', 2);
plot(theta_i_range, db(Pr_Theta_beta2_matrix(:, idx_best_RMSE), 'power'), 'r-s', 'LineWidth', 2);
xlabel('Incidence Angle (°)', 'FontSize', 12);
ylabel('Power (dB)', 'FontSize', 12);
title(sprintf('p_r and P_{IT} vs. Incidence angle at Optimal \\beta_2 = %.2f (Rs = %.1f m)', ...
    beta2_vals(idx_best_RMSE), Rs));
legend({'P_{IT}', 'p_r'}, 'Location', 'best');
set(gca, 'FontSize', 10);

%% ABSe

num_theta_ABSe = size(ABSe_cube, 1); 
best_beta2_idx_ABSe = zeros(num_theta_ABSe, 1);
best_Pr_ABSe = zeros(num_theta_ABSe, 1);
best_PIT_ABSe = zeros(num_theta_ABSe, 1);

for Theta_idx_ABSe = 1:num_theta_ABSe
    [~, best_idx_ABSe] = min(ABSe_cube(Theta_idx_ABSe,:,1)); % index of min ABSe
    best_beta2_idx_ABSe(Theta_idx_ABSe) = best_idx_ABSe;
    best_Pr_ABSe(Theta_idx_ABSe) = Pr_Theta_beta2_matrix(Theta_idx_ABSe, best_idx_ABSe);
    best_PIT_ABSe(Theta_idx_ABSe) = P_IT_matrix_all(Theta_idx_ABSe, best_idx_ABSe);
end

figure;
plot(theta_i_range, db(best_Pr_ABSe, 'power'), '-o', 'LineWidth', 2, 'MarkerSize', 3);
hold on;
plot(theta_i_range, db(best_PIT_ABSe, 'power'), '-s', 'LineWidth', 2, 'MarkerSize', 3);
legend({'Pr', 'P\_IT'}, 'Location', 'best');
xlabel('Incidence Angle (°)');
ylabel('Power (dB)');
title('Best \beta_2 (ABSe) for Each \theta');
grid on;

% Add labels
for i = 1:length(theta_i_range)
    beta2_val_text = num2str(beta2_vals(best_beta2_idx_ABSe(i)), '%.2f');
    text(theta_i_range(i), db(best_Pr_ABSe(i), 'power') + 0.5, beta2_val_text, ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k'); 
end
%% ABSe

end
end


%% RADAR Equation
function [Pr, coeff_matrix] = radarEquation(R0, Rs, lambda, sigma_z, beta_c, ...
                                   Pt, Ny_surf, Nx_surf, ...
                                   freq, SMvol, sigma_z_cm, R0_slant, Rs_slant, ...
                                   beta0, beta2,theta_i, phi_i, theta_s_vector, phi_s_vector,surfParam,TRAp,sigmaParam_Betac,geoSYSp,Rx_Gain_All,RxPmod,Cords,Area)

R2 = RxPmod;
R2_vector = R2(:);

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
theta_s_vector_filtered(theta_s_vector > 10000) = NaN;

theta_s_rad = deg2rad(theta_s_vector_filtered);
phi_s_rad = deg2rad(phi_s_vector);

coeff_vector = NaN(size(theta_s_rad));

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
R1 = Area.R1;

dA        = (Area.x_vec(2)-Area.x_vec(1)) * (Area.y_vec(2)-Area.y_vec(1));
integrand = HSAVERSTXgain(TRAp,Rs,Cords) .* Rx_Gain_All.RxGain ...
        .* coeff_matrix ./ (R1.^2 .* R2.^2);

Pr = ((Pt * lambda^2)/(4*pi)^3) * sum(integrand(:),'omitnan') * abs(dA);
end

%% HSAVERS TX gain
function TxGain = HSAVERSTXgain(TRAp,Rs,Cords)
GPS_theta         = TRAp.theta_Txpattern;
GPS_phi           = TRAp.phi_Txpattern; 
GPSPattern_linear = TRAp.GainTx_linear(:,:,1);
TRAp.TxAttitude_angles(3) = TRAp.TxAttitude_angles(3) +  TRAp.Yaw_TX;
[Txmat_ECEF2OF, Txmat_OF2BF, Txmat_BFtoAF] = matrixFrameGen(TRAp.TxECEF,  TRAp.VTxECEF, TRAp.Tx_nomAtt_angles, TRAp.TxAttitude_angles, TRAp.TxEuler_angles);

RxECEF_Rs = Cords.RxECEF_Rs;
TRx_K = (RxECEF_Rs - TRAp.TxECEF);
TRx_OF = Txmat_ECEF2OF*TRx_K;
TRx_BF = Txmat_OF2BF*TRx_OF;
TRx_AF = Txmat_BFtoAF*TRx_BF;

%piercing angles
thetaGPS_pier = acosd(dot(TRx_AF./norm(TRx_AF), [0,0,1]));
phiGPS_pier   = atan2d(TRx_AF(2), TRx_AF(1));
[~,ind_phiTransm] = min(abs(phiGPS_pier - GPS_phi));
[~,ind_thTransm]  = min(abs(thetaGPS_pier    - GPS_theta));
TxGain = GPSPattern_linear(ind_phiTransm, ind_thTransm);

end


%% RX gain
function [RxGain,Gain_LHCP_onsurf_lin_atSP,Gain_RHCP_onsurf_lin_atSP] = HSAVERSRXgain(TRAp,Xplot,Yplot,Nx_surf,Ny_surf,Rs,Cords,Area,geoSYSp)

pSGR(:,:,1) = Area.pSGR(:,:,1);
pSGR(:,:,2) = Area.pSGR(:,:,2);
pSGR(:,:,3) = Area.pSGR(:,:,3);     
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

if geoSYSp.Rx_Attitude_Variation_RPY(1) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(2) == 1 || geoSYSp.Rx_Attitude_Variation_RPY(3) == 1
RxAttitude = [TRAp.Rx_roll,TRAp.Rx_pitch,0];
else
RxAttitude =  TRAp.Rx_nomAtt_angles;
end

[Amat_ECEF2OF, Amat_OF2BF, Amat_BFtoAF] = matrixFrameGen(RxECEF, vRxECEF, RxAttitude, TRAp.RxAttitude_angles, TRAp.RxEuler_angles);

        
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

function [x_vec, y_vec, X, Y, Nx_surf, Ny_surf, Area] = precomputeGrid(geoSYSp, Rs, TRAp_Betac,theta_i)

% theta_i = TRAp_Betac.incAngle_enuRot;
Xmax = round(sqrt(2 * 2.99792458e2 * geoSYSp.Dly_RNGmax * Rs * cosd(theta_i)) ...
        / cosd(theta_i)^2 / 1000) * 1000;

% if Xmax > 10000 
% dx = 500;
% end 

% Generate X and Y vectors
XplotPos = 0.0:geoSYSp.dx:Xmax;
XplotNeg = -(geoSYSp.dx:geoSYSp.dx:Xmax);
Xplot = [fliplr(XplotNeg), XplotPos];
x_vec = Xplot;
Nx_surf = length(Xplot);

YplotPos = 0.0:geoSYSp.dy:Xmax;
YplotNeg = -YplotPos(2:end);
Yplot = [fliplr(YplotNeg), YplotPos]';
y_vec = Yplot;
Ny_surf = length(Yplot);

% Create 2D meshgrids
[X, Y] = meshgrid(x_vec, y_vec);
[xp, yp] = meshgrid(x_vec, y_vec);
zp = zeros(size(xp));

% Create 3D surface grid
pSGR = cat(3, xp, yp, zp);

% Compute R1 (slant range from Tx to surface)
[~, ~, R1_scalar] = xyz2llh(TRAp_Betac.TxECEF(1), TRAp_Betac.TxECEF(2), TRAp_Betac.TxECEF(3));
R1 = repmat(R1_scalar, Ny_surf, Nx_surf, 1);

% Store everything in Area struct
Area = struct( ...
    'Xmax_Range', Xmax, ...
    'Xmin_Range', -Xmax, ...
    'Xplot', Xplot, ...
    'Yplot', Yplot, ...
    'x_vec', x_vec, ...
    'y_vec', y_vec, ...
    'Nx_surf', Nx_surf, ...
    'Ny_surf', Ny_surf, ...
    'X', X, ...
    'Y', Y, ...
    'xp', xp, ...
    'yp', yp, ...
    'zp', zp, ...
    'pSGR', pSGR, ...
    'R1', R1 ...
);
end




% function [Cords, pSGR] = computeCords(Rs, TRAp_Betac, geoSYSp, Area, Rs_range,needed_inc)
% 
%  % Precompute grids
% regGridSim_z = zeros(Area.Ny_surf, Area.Nx_surf);
% xp = Area.X;
% yp = Area.Y;
% zp = regGridSim_z;
% pSGR = cat(3, xp, yp, zp);
% 
% % Convert Rx to ECEF for this Rs
% if numel(Rs_range) == 1
% [lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, ~] = xyz2llh(TRAp_Betac.RxECEF(1), TRAp_Betac.RxECEF(2), TRAp_Betac.RxECEF(3));
% lat_Rx_Temp_Rs = rad2deg(lat_Rx_Temp_Rs);
% lon_Rx_Temp_Rs = rad2deg(lon_Rx_Temp_Rs);
% RxLLA_Rs_temp = [lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, Rs];
% RxECEF_Rs = lla2ecef(RxLLA_Rs_temp)';
% else
% RxECEF_Rs = TRAp_Betac.RxECEF;
% end
% 
% %% find fake GPS position to have desired incidence angle,
% TxECEF_orig = TRAp_Betac.TxECEF; % Original TxECEF
% best_TxECEF = TxECEF_orig;
% best_error = inf;
% 
% [Tx_lat_orig, Tx_lon_orig, Tx_alt_orig] = xyz2llh(TRAp_Betac.TxECEF(1),TRAp_Betac.TxECEF(2),TRAp_Betac.TxECEF(3));
% Tx_lat_orig = rad2deg(Tx_lat_orig);
% Tx_lon_orig = rad2deg(Tx_lon_orig);
% 
% %%
% objective_xyz = @(vars) abs(calcIncAngle(vars(:), RxECEF_Rs) - needed_inc);
% 
% tolerances = [1e-3, 1e-2, 1e-1, 1];
% for tol = tolerances
%     options = optimset('TolX', tol, 'Display', 'off', 'MaxIter', 200);
%     [candidate, fval] = fminsearch(objective_xyz, best_TxLLH, options);
%     if fval < best_error
%         best_error = fval;
%         best_TxLLH = candidate(:);
%     else
%     end
% end
% 
% TxLLH_new = best_TxLLH;
% 
% %%
% [TxECEF_new(1),TxECEF_new(2),TxECEF_new(3)] = llh2xyz(TxLLH_new(1),TxLLH_new(2),TxLLH_new(3));
% % RECOMPUTE COORDINATES WITH NEW TxECEF ---
% OENUinLLA = Specular_Point_Position_Calculation(RxECEF_Rs', TxECEF_new');
% OENUinECEF = lla2ecef(OENUinLLA)';
% 
% TxENU_Rs = ecef2enu_dc(TxECEF_new, OENUinECEF);
% RxENU_Rs = ecef2enu_dc(RxECEF_Rs, OENUinECEF);
% 
% phiTx = atan2(TxENU_Rs(2), TxENU_Rs(1));
% phi_ENU2local = pi - phiTx;
% Txlocal_rs = rotation_clockwise(TxENU_Rs, phi_ENU2local, 'rad');
% Rxlocal_rs = rotation_clockwise(RxENU_Rs, phi_ENU2local, 'rad');
% 
% %% double check
% 
% gamma_new = -atan(Txlocal_rs(3) / Txlocal_rs(1));
% incAngle_new = 90 - rad2deg(gamma_new);
% 
% %%
% 
% % Output structure
% Cords.Txlocal_rs     = Txlocal_rs;
% Cords.Rxlocal_rs     = Rxlocal_rs;
% Cords.OENUinECEF     = OENUinECEF;
% Cords.phi_ENU2local  = phi_ENU2local;
% Cords.TxENU_Rs       = TxENU_Rs;
% Cords.RxENU_Rs       = RxENU_Rs;
% Cords.RxECEF_Rs      = RxECEF_Rs;
% Cords.calculated_inc = incAngle_new;
% Cords.Needed_inc = needed_inc;
% 
% end
% 
% % --- HELPER FUNCTION ---
% function incAngle = calcIncAngle(TxECEF, RxECEF)
%     OENUinLLA = Specular_Point_Position_Calculation(RxECEF', TxECEF');
%     OENUinECEF = lla2ecef(OENUinLLA)';
% 
%     TxENU = ecef2enu_dc(TxECEF, OENUinECEF);
%     phiTx = atan2(TxENU(2), TxENU(1));
%     phi_ENU2local = pi - phiTx;
%     Txlocal = rotation_clockwise(TxENU, phi_ENU2local, 'rad');
% 
%     gamma = -atan(Txlocal(3) / Txlocal(1));
%     incAngle = 90 - rad2deg(gamma);
% end
% 

function [Cords, pSGR] = computeCords(Rs, TRAp_Betac, geoSYSp, Area, Rs_range, needed_inc)

    % Precompute grids
    regGridSim_z = zeros(Area.Ny_surf, Area.Nx_surf);
    xp = Area.X;
    yp = Area.Y;
    zp = regGridSim_z;
    pSGR = cat(3, xp, yp, zp);

    % Convert Rx to ECEF for this Rs
    if numel(Rs_range) == 1
        [lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, ~] = xyz2llh( ...
            TRAp_Betac.RxECEF(1), TRAp_Betac.RxECEF(2), TRAp_Betac.RxECEF(3));
        lat_Rx_Temp_Rs = rad2deg(lat_Rx_Temp_Rs);
        lon_Rx_Temp_Rs = rad2deg(lon_Rx_Temp_Rs);
        RxLLA_Rs_temp = [lat_Rx_Temp_Rs, lon_Rx_Temp_Rs, Rs];
        RxECEF_Rs = lla2ecef(RxLLA_Rs_temp)';
    else
        RxECEF_Rs = TRAp_Betac.RxECEF;
    end

    [Tx_lat_orig, Tx_lon_orig, Tx_alt_orig] = xyz2llh( ...
        TRAp_Betac.TxECEF(1), TRAp_Betac.TxECEF(2), TRAp_Betac.TxECEF(3));
    Tx_lat_orig = rad2deg(Tx_lat_orig);
    Tx_lon_orig = rad2deg(Tx_lon_orig);

    TxLLH_orig = [Tx_lat_orig, Tx_lon_orig, Tx_alt_orig];

    %% Define objective function for lat/lon (altitude fixed)
    objective_llh = @(vars) calcIncAngle_LLh([vars(:).' TxLLH_orig(3)], RxECEF_Rs, needed_inc);

    lat_lb = -90;   lat_ub = 90;
    lon_lb = -180;  lon_ub = 180;

    %% fminsearchbnd
    options = optimset('Display', 'off', 'TolX', 1e-3, 'MaxIter', 300);
    [candidate, fval] = fminsearchbnd(objective_llh, TxLLH_orig(1:2), ...
                                      [lat_lb, lon_lb], [lat_ub, lon_ub], options);

    best_TxLLH = [candidate(:).' TxLLH_orig(3)];

    %% Convert LLH back to ECEF
    [TxECEF_new(1),TxECEF_new(2),TxECEF_new(3)] = llh2xyz(best_TxLLH(1), best_TxLLH(2), best_TxLLH(3));

    %%Recompute coordinates with new TxECEF
    OENUinLLA = Specular_Point_Position_Calculation(RxECEF_Rs(:).', TxECEF_new(:).');
    OENUinECEF = lla2ecef(OENUinLLA)';

    TxENU_Rs = ecef2enu_dc(TxECEF_new, OENUinECEF);
    RxENU_Rs = ecef2enu_dc(RxECEF_Rs, OENUinECEF);

    phiTx = atan2(TxENU_Rs(2), TxENU_Rs(1));
    phi_ENU2local = pi - phiTx;
    Txlocal_rs = rotation_clockwise(TxENU_Rs, phi_ENU2local, 'rad');
    Rxlocal_rs = rotation_clockwise(RxENU_Rs, phi_ENU2local, 'rad');

    %% Double check incidence angle
    gamma_new = -atan(Txlocal_rs(3) / Txlocal_rs(1));
    incAngle_new = 90 - rad2deg(gamma_new);

    %% Output structure
    Cords.Txlocal_rs     = Txlocal_rs;
    Cords.Rxlocal_rs     = Rxlocal_rs;
    Cords.OENUinECEF     = OENUinECEF;
    Cords.phi_ENU2local  = phi_ENU2local;
    Cords.TxENU_Rs       = TxENU_Rs;
    Cords.RxENU_Rs       = RxENU_Rs;
    Cords.RxECEF_Rs      = RxECEF_Rs;
    Cords.calculated_inc = incAngle_new;
    Cords.Needed_inc     = needed_inc;
    Cords.TxLLH_new      = best_TxLLH;

end

function error_val = calcIncAngle_LLh(TxLLH, RxECEF_Rs, needed_inc)
    [TxECEF(1),TxECEF(2),TxECEF(3)] = llh2xyz(TxLLH(1), TxLLH(2), TxLLH(3));
    incAngle = calcIncAngle(TxECEF, RxECEF_Rs);
    error_val = abs(incAngle - needed_inc);
end

function incAngle = calcIncAngle(TxECEF, RxECEF)
    OENUinLLA = Specular_Point_Position_Calculation(RxECEF(:).', TxECEF(:).');
    OENUinECEF = lla2ecef(OENUinLLA)';
    TxENU = ecef2enu_dc(TxECEF, OENUinECEF);
    phiTx = atan2(TxENU(2), TxENU(1));
    phi_ENU2local = pi - phiTx;
    Txlocal = rotation_clockwise(TxENU, phi_ENU2local, 'rad');
    gamma = -atan(Txlocal(3) / Txlocal(1));
    incAngle = 90 - rad2deg(gamma);
end
