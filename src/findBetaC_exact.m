function beta_c = findBetaC_exact(R0, Rs, lambda, Pt_dBW, sigma_z_cm, ...
                                  ell_x_cm, ell_y_cm, Gain2D, Ny, Nx)
% FINDBETAC_EXACT  Exact implementation of Eqs.(3),(10),(13) from Comite et al.
%
% beta_c = findBetaC_exact(R0, Rs, lambda, Pt_dBW, sigma_z_cm, ...
%                          ell_x_cm, ell_y_cm, Gain2D, Ny, Nx)
%
% Inputs:
%   R0, Rs         – Tx/Rx heights above surface [m]
%   lambda         – wavelength [m]
%   Pt_dBW         – transmit EIRP [dBW]
%   sigma_z_cm     – RMS surface roughness [cm]
%   ell_x_cm       – surface corr-length in x [cm]
%   ell_y_cm       – surface corr-length in y [cm]
%   Gain2D         – nAz×nInc measured LHCP gains (linear)
%   Ny, Nx         – patch resolution in y/x
%
% Output:
%   beta_c         – optimal coherent beamwidth [deg]

  %% 0) Convert units
  Pt      = 10.^(Pt_dBW/10);    % dBW→W
  sigma_z = sigma_z_cm/100;     % cm→m
  ell_x   = ell_x_cm  /100;     % cm→m
  ell_y   = ell_y_cm  /100;     % cm→m

  %% 1) Build angle axes for Gain2D (assumes full 0–360,0–90 coverage)
  [nAz, nInc] = size(Gain2D);
  azi_rad     = linspace(0,360,nAz)*pi/180;  % 1×nAz
  inc_rad     = linspace(0, 90,nInc)*pi/180;  % 1×nInc

  %% 2) Create true 2D cubic interpolant
  F2D = griddedInterpolant({azi_rad,inc_rad}, Gain2D, 'cubic', 'none');
  Dm  = @(th,ph) max(0, F2D(ph,th));

  %% 3) Define specular angles θ_s [0°→60°]
  theta_s = linspace(0,60,61)*pi/180;  % rad
  N       = numel(theta_s);

  %% 4) Compute image‐theory reference P_IT (Eq.3)
  P_IT = zeros(1,N);
  for i = 1:N
    th    = theta_s(i);
    a_pp  = reflectionCoeff(th);
    P_IT(i) = Pt * Dm(th,0)^2 ...
              * lambda^2 * abs(a_pp)^2 ...
              / ((4*pi)^2 * (R0+Rs)^2);
  end

  %% 5) Precompute RMSE curve for diagnostics
  beta_test = linspace(0.01,5,200);
  rmse_vals = arrayfun(@(b) computeRMSE(b), beta_test);

  figure;
  plot(beta_test, rmse_vals, 'LineWidth',1.5);
  grid on;
  xlabel('\beta_c [deg]'); ylabel('RMSE');
  title('RMSE vs. \beta_c (Exact Eq.10)');

  %% 6) Optimize βc over a bracket likely containing the minimum
  % find approximate minimum location from the curve
  [~,idx_min] = min(rmse_vals);
  b0 = max(0.01, beta_test(idx_min)-1);
  b1 = min(5, beta_test(idx_min)+1);

  opts = optimset('TolX',1e-4,'Display','off');
  beta_c = fminbnd(@computeRMSE, b0, b1, opts);

  fprintf('Optimal β_c = %.4f° (search window [%.2f,%.2f])\n', beta_c, b0, b1);


  %% Nested RMSE function
  function err = computeRMSE(beta_c)
    Pr = arrayfun(@(th) Pradar(th, beta_c), theta_s);
    err = sqrt( mean( (Pr - P_IT).^2 ) );
  end

  %% Nested radar‐equation implementing full Eq.(10)&(13)
  function Pr = Pradar(th_s, beta_c)
    k  = 2*pi/lambda;
    kz = k*sin(th_s);

    % Patch coordinates
    L   = sqrt(lambda * Rs);
    x   = linspace(-L, L, Nx);
    y   = linspace(-L, L, Ny);
    [X,Y] = meshgrid(x,y);

    % Ranges
    R1 = sqrt(R0^2 + X.^2 + Y.^2);
    R2 = sqrt(Rs^2 + X.^2 + Y.^2);

    % Local angles
    th1 = atan2(sqrt(X.^2+Y.^2), R0);
    th2 = atan2(sqrt(X.^2+Y.^2), Rs);
    ph1 = atan2(Y, X);
    ph2 = ph1;  % specular reciprocity

    % Fresnel at specular point
    a_pp = reflectionCoeff(th_s);

    % Full‐form Eq.(10) parameters
    dx    = th1 - th_s;
    dy    = th2 - th_s;
    kx    = k*sin(th1);
    ky    = k*sin(th2);
    delta = (ell_x^2 * kx + ell_y^2 * ky)/2;

    a3 = 1/(R0*(beta_c*pi/180)) + (k^2 * ell_y^2)/4;
    a4 = 1/(Rs*(beta_c*pi/180)) + (k^2 * ell_x^2)/4;

    sigma_c = (k^2 * abs(a_pp)^2) ./ (4*sqrt(a3.*a4)) ...
            .* exp(-kz.^2 * sigma_z^2) ...
            .* exp(-k^2*dy.^2./(4*a3)) ...
            .* exp(-k^2*dx.^2./(4*a4)) ...
            .* exp(-(kx.*ky).*delta./(2*a3.*a4));

    % Discretized integral (Eq.13)
    dA        = (x(2)-x(1)) * (y(2)-y(1));
    integrand = Dm(th1,ph1).*Dm(th2,ph2) ...
              .* sigma_c ./ (R1.^2 .* R2.^2);

    Pr = Pt * (lambda^2/(4*pi)^3) * sum(integrand(:)) * dA;
  end
end

%% Eq.2 Fresnel reflection coefficient
function a_pp = reflectionCoeff(theta)
  eps_r = 15 - 1i*0;
  ct    = cos(theta);
  s     = sqrt(eps_r - sin(theta).^2);
  a_pp  = (eps_r*ct - s) ./ (eps_r*ct + s);
end
