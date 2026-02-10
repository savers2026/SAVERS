function beta_c_opt = findBetaC(R0, Rs, lambda, Pt, sigma_z, ...
                                 TRApGain_rx, Ny_surf, Nx_surf, TRApGain_tx)
% findBetaC Compute optimal effective beamwidth β_c (Simplified)

% Inputs:
%   R0         -- transmitter height [m]
%   Rs         -- receiver height [m]
%   lambda     -- wavelength [m]
%   Pt         -- Tx power [W]
%   sigma_z    -- RMS surface roughness [m]
%   TRApGain_rx   -- RX Antenna gain pattern
%   Ny_surf
%   Nx_surf
%   TRApGain_tx  -- TX Antenna gain pattern
%
% Output:
%   beta_c_opt -- optimal β_c [deg] minimizing RMSE over θ_s∈[0,60°]

Pt = db2pow(Pt); %W
sigma_z = sigma_z/100; %m

  inc_rad_rx = linspace(0,90,size(TRApGain_rx,2)) * pi/180;
  azi_rad_rx = linspace(0,360,size(TRApGain_rx,1)) * pi/180;
  F2D_rx = griddedInterpolant({azi_rad_rx, inc_rad_rx}, TRApGain_rx, ...
                           'cubic','none');
  inc_rad_tx = linspace(0,90,size(TRApGain_tx,2)) * pi/180;
  azi_rad_tx = linspace(0,360,size(TRApGain_tx,1)) * pi/180;
  F2D_tx = griddedInterpolant({azi_rad_tx, inc_rad_tx}, TRApGain_tx, ...
                           'cubic','none');

  D_rx = @(theta,phi) max(0, F2D_rx(phi,theta));
  D_tx = @(theta,phi) max(0, F2D_tx(phi,theta));

  theta_s = linspace(0,60,61) * pi/180;
  N       = numel(theta_s);

  %Compute image theory power P_IT(θ_s)
  P_IT = zeros(1,N);
  for i = 1:N
    th       = theta_s(i);
    a_pp     = reflectionCoeff(th);
    P_IT(i)  = Pt * D_rx(th,0) * D_tx(th,0) ...         % assume azimuth=0 (symmetric)
               * lambda^2 * abs(a_pp)^2 ...
               / ((4*pi)^2 * (R0 + Rs)^2);
  end

  %RMSE over β_c
  rmseFun = @(beta_c) sqrt( mean( ...
      arrayfun(@(th, Pit) ...
        (radarEquation(R0,Rs,lambda,sigma_z,beta_c,th,Pt,D_rx,D_tx,Ny_surf,Nx_surf) ...
         - Pit).^2, ...
        theta_s, P_IT) ) );

  %Minimize RMSE over β_c
  opts       = optimset('TolX',1e-4,'Display','off');
  beta_c_opt = fminbnd(rmseFun, 0.0001, 100, opts);

end

%% --- Radar Equatuin ---------------------------------------------

function Pr = radarEquation(R0, Rs, lambda, sigma_z, beta_c, ...
                            theta_s, Pt, D_rx, D_tx, Ny_surf, Nx_surf)
  % Integrate coherent NRCS over the surface patch

  % Patch around specular point of size Ny_surf×Nx_surf
  L      = sqrt(lambda * Rs);
  x_vec  = linspace(-L, L, Nx_surf);
  y_vec  = linspace(-L, L, Ny_surf);
  [X,Y]  = meshgrid(x_vec, y_vec);

  % Ranges to each element
  R1 = sqrt(R0^2 + X.^2 + Y.^2);
  R2 = sqrt(Rs^2 + X.^2 + Y.^2);

  % Local angles
  theta1 = atan2(sqrt(X.^2 + Y.^2), R0);
  theta2 = atan2(sqrt(X.^2 + Y.^2), Rs);
  phi1   = atan2(Y, X);
  phi2   = phi1;

  % Coherent NRCS parameters from Eq. (10)
  g0     = 1 / (R0 * (beta_c*pi/180));
  gs     = 1 / (Rs * (beta_c*pi/180));
  kz     = (4*pi/lambda) * sin(theta_s);
  delta2 = (theta1 - theta_s).^2 + (theta2 - theta_s).^2;

  sigma_c = abs(reflectionCoeff(theta_s)).^2 ...
            * (1/(g0*gs)) ...
            * exp(-kz.^2 * sigma_z^2) ...
            * exp(-delta2 * (g0 + gs));

  dA        = (x_vec(2)-x_vec(1)) * (y_vec(2)-y_vec(1));
  integrand = D_tx(theta1,phi1) .* D_rx(theta2,phi2) ...
            .* sigma_c ./ (R1.^2 .* R2.^2); % .* sigma_c

  Pr = Pt * (lambda^2/(4*pi)^3) * sum(integrand(:)) * dA;
end

function a_pp = reflectionCoeff(theta)
  % Fresnel reflection coefficient for flat surface
  eps_r    = 13;
  cos_t    = cos(theta);
  sqrtTerm = sqrt(eps_r - sin(theta).^2);
  a_pp     = (eps_r*cos_t - sqrtTerm) ./ (eps_r*cos_t + sqrtTerm);
end

