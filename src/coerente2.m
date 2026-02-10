function coeff = coerente(f, r, sstd, R0, R2, beta0rad, beta2rad, t, p, ts, ps)
% coerente: Compute the coherent part of the bistatic scattering coefficient
% Inputs:
%   f         - frequency in Hz
%   r         - complex reflection coefficient
%   sstd      - RMS height of Gaussian rough surface (m)
%   R0        - slant range of transmit antenna (m)
%   R2        - slant range of receive antenna (m)
%   beta0rad  - transmit antenna beamwidth (rad)
%   beta2rad  - receive antenna beamwidth (rad)
%   t, p      - incidence angles theta and phi (rad)
%   ts, ps    - scattering angles theta_s and phi_s (rad)
% Output:
%   coeff     - coherent part of bistatic scattering coefficient


c   = 3e8;             
lam = c / f;           
k   = 2 * pi / lam;    
j   = 1i;

g0 = 1 / (R0 * beta0rad);
gs = 1 / (R2 * beta2rad);

theta   = t;
phi     = p;
theta_s = ts;
phi_s   = ps;

kx = k * (sin(theta_s)*cos(phi_s) - sin(theta));
if phi_s == pi
    ky = 0;
else
    ky = k * sin(theta_s) * sin(phi_s);
end
kz = k * (cos(theta) + cos(theta_s));

ax0  = cos(theta)^2;
ax2  = 1 - (sin(theta_s)^2)*(cos(phi_s)^2);
ay2  = 1 - (sin(theta_s)^2)*(sin(phi_s)^2);
axy2 = 2 * sin(phi_s) * cos(phi_s) * sin(theta_s)^2;

alphas = (cos(theta_s)^2)*(cos(phi_s)^2) + sin(phi_s)^2;
betas  = (cos(theta_s)^2)*(sin(phi_s)^2) + cos(phi_s)^2;
gammas = -2 * sin(phi_s)*cos(phi_s)*sin(theta_s)^2;

alpha = (g0^2)*cos(theta)^2 + (gs^2)*alphas;
beta  = g0^2 + (gs^2)*betas;
gamma = (gs^2)*gammas;

a1 = 2 * beta;
a2 = 2 * alpha - (gamma^2)/a1;

a3 = beta/2 + ((k^2)/(4*a1) + (k^2)*(gamma^2)/(4*a2*a1^2)) * (1/(R0^2) + (ay2^2)/(R2^2) + 2*ay2/(R0*R2)) + ...
     (k^2)/(4*a2) * (axy2^2)/(4*R2^2) + ((k^2)*gamma/(2*a2*a1)) * (axy2/(2*R0*R2) + axy2*ay2/(2*R2^2));

delta3 = gamma/2 - ((k^2)/(4*a1) + (k^2)*(gamma^2)/(4*a2*a1^2)) * (axy2/(R0*R2) + ay2*axy2/(R2^2)) - ...
         (k^2)/(4*a2) * (ax0*axy2/(R0*R2) + ax2*axy2/(R2^2)) - ((k^2)*gamma/(2*a2*a1)) * ...
         (ax0/(R0^2) + ax2/(R0*R2) + ax0*ay2/(R0*R2) + ax2*ay2/(R2^2) + (axy2^2)/(4*R2^2));

a4 = alpha/2 + ((k^2)/(4*a1) + (k^2)*(gamma^2)/(4*a2*a1^2)) * (axy2^2)/(4*R2^2) + ...
     (k^2)/(4*a2) * ((ax0^2)/(R0^2) + (ax2^2)/(R2^2) + 2*ax0*ax2/(R0*R2)) + ...
     ((k^2)*gamma/(2*a2*a1)) * (ax0*axy2/(2*R0*R2) + ax2*axy2/(2*R2^2)) - delta3^2/(4*a3);

b4 = -j * kx + j * ky * delta3/(2*a3);

a0 = abs(r) * (cos(theta) + cos(theta_s));

coeff = (k^2)*(a0^2)/4 * sqrt(1/a3) * sqrt(1/a4) * exp(-(kz^2)*sstd^2) * ...
        exp(-0.25*(ky^2)/a3) * exp((b4^2)/(4*a4));
end
