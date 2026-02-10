function [PowerDD, WAF] = MapinDDandWAF(ACF, Refl_Surf, DDM, Tc)

%% ==============================================================================================
% To generate the power of the reflected signal in Delay-Doppler domain from the surface domain
% INPUT:
  % ACF: structure including the ACF parameters
  % Refl_Surf: structure including the surface parameters:
    % - delay_surface_grid: delay of each surface cell
    % - Doppler_surface_grid: Doppler of each surface cell
    % - power_surface_grid: power of reflected signal of each surface cell
  % DDM: structure including the DDM parameters, e.g.
    % Delay_resolution, DDM_delay, Doppler_resolution, DDM_Doppler
  % Tc: coherent integration time
% Output:
  % WAF: Woodward's ambiguity function adapted to the DDM resolution
  % PowerDD: reflected signal power in Delay and Doppler Domain
%=========================================================================================================

%% ACF and Reflection Surface Paramters
ACF_range_delay = ACF.ACF_range_delay;
ACF = ACF.ACF;
delay_surface_grid = Refl_Surf.delay_surface_grid;
Doppler_surface_grid = Refl_Surf.Doppler_surface_grid;
power_surface_grid = Refl_Surf.power_surface_grid;
Delay_resolution = DDM.Delay_resolution;
DDM_delay = DDM.DDM_delay;
Doppler_resolution = DDM.Doppler_resolution;
DDM_Doppler = DDM.DDM_Doppler;


%% delay and Doppler range and resolution with interpolation factors
delay_itp_factor = 1.0;           % interpolation factor of the delay bins when generating the power in delay Doppler domain (before the convolution with WAF)
Doppler_itp_factor = 5.0;         % interpolation factor of the Doppler bins when generating the power in delay Doppler domain (before the convolution with WAF)
chiplength = 299792458.0/1.023e6;

samples_per_chip = ceil(chiplength/Delay_resolution*delay_itp_factor);
min_delay_bin = min(DDM_delay/Delay_resolution);
max_delay_bin = max(DDM_delay/Delay_resolution);
delay_resolution_itp = Delay_resolution/delay_itp_factor;
delay_Sigma = ((min_delay_bin*delay_itp_factor - 2.0*samples_per_chip):(max_delay_bin*delay_itp_factor +  2.0*samples_per_chip))*delay_resolution_itp;

min_Doppler_bin = min(DDM_Doppler/Doppler_resolution);
max_Doppler_bin = max(DDM_Doppler/Doppler_resolution);
Doppler_resolution_itp = Doppler_resolution/Doppler_itp_factor;

min_Doppler = (min_Doppler_bin*Doppler_itp_factor - 3.0/Tc/Doppler_resolution_itp)*Doppler_resolution_itp;
max_Doppler = (max_Doppler_bin*Doppler_itp_factor + 3.0/Tc/Doppler_resolution_itp)*Doppler_resolution_itp;

tmp = min_Doppler*Tc;
if tmp < 0.0
    tmp = fix(tmp) - 1;
else
    tmp = fix(tmp) + 1;
end
min_Doppler = tmp/Tc;

tmp = max_Doppler*Tc;
if tmp < 0.0
    tmp = fix(tmp) - 1;
else
    tmp = fix(tmp) + 1;
end
max_Doppler = tmp/Tc;
Doppler_Sigma = min_Doppler:Doppler_resolution_itp:max_Doppler;

[c, sp_idx_delay] = min(abs(delay_Sigma));
[c, sp_idx_Doppler] = min(abs(Doppler_Sigma));
DDM_idx_delay = (sp_idx_delay + min_delay_bin*delay_itp_factor):delay_itp_factor:(sp_idx_delay + max_delay_bin*delay_itp_factor);
DDM_idx_Doppler = (sp_idx_Doppler + min_Doppler_bin*Doppler_itp_factor):Doppler_itp_factor:(sp_idx_Doppler + max_Doppler_bin*Doppler_itp_factor);


%%
% interpolations of the delay, Doppler and power over each surface grid
% [s1, s2] = size(delay_surface_grid);
% [X,Y] = meshgrid(1:s1);
% [Xq,Yq] = meshgrid(1:0.1:s1);
% delay_surface_grid_itp = interp2(X,Y,delay_surface_grid,Xq,Yq, 'cubic');
% Doppler_surface_grid_itp = interp2(X,Y,Doppler_surface_grid,Xq,Yq, 'cubic');
% power_surface_grid_itp = interp2(X,Y,power_surface_grid,Xq,Yq, 'cubic');
% power_surface_grid_itp = power_surface_grid_itp*0.1*0.1;
% Doppler_surface_grid = Doppler_surface_grid_itp;
% power_surface_grid = power_surface_grid_itp;
% delay_surface_grid = delay_surface_grid_itp;


%% Generate the reflected power in Delay Doppler domain
power_in_DD = zeros(length(Doppler_Sigma), length(delay_Sigma));

for i = 1:length(Doppler_Sigma)
    valid_Doppler_grid = (Doppler_surface_grid >= (Doppler_Sigma(i) - Doppler_resolution_itp/2.0)) & (Doppler_surface_grid < (Doppler_Sigma(i) + Doppler_resolution_itp/2.0));
    delay_valid = delay_surface_grid(valid_Doppler_grid);
    power_valid = power_surface_grid(valid_Doppler_grid);
    for j = 1:length(delay_Sigma)
        valid_DD_grid = (delay_valid >= (delay_Sigma(j) - delay_resolution_itp/2.0)) & (delay_valid < (delay_Sigma(j) + delay_resolution_itp/2.0));
        power_in_DD(i,j) = sum(power_valid(valid_DD_grid));
    end
end

PowerDD = {};
PowerDD.power_in_DD = power_in_DD;
PowerDD.delay_resolution_itp = delay_resolution_itp;
PowerDD.Doppler_resolution_itp = Doppler_resolution_itp;
PowerDD.delay_Sigma = delay_Sigma;
PowerDD.Doppler_Sigma = Doppler_Sigma;
PowerDD.DDM_idx_delay = DDM_idx_delay;
PowerDD.DDM_idx_Doppler = DDM_idx_Doppler;


%% WAF with the delay and Doppler resolution compatible with power_in_DD
n_sample = ceil(2.5*chiplength/delay_resolution_itp);
Delay_WAF = (-n_sample:n_sample)*delay_resolution_itp;
Doppler_WAF = -4.0/Tc:Doppler_resolution_itp:4.0/Tc;
WAF0 = interp1(ACF_range_delay, ACF, Delay_WAF, 'PCHIP');
WAF1 = sinc(Tc*Doppler_WAF).^2.0;

Delta_Delay = (0:n_sample)*Delay_resolution;
ACFbyACF = zeros(length(Delta_Delay), length(Delay_WAF));
for i = 1:length(Delta_Delay)
    d_delay = Delta_Delay(i);
    ACF1 = interp1(ACF_range_delay, ACF, Delay_WAF, 'PCHIP', 0);
    ACF2 = interp1(ACF_range_delay, ACF, Delay_WAF + d_delay, 'PCHIP', 0);
    ACFbyACF(i,:) = ACF1.*ACF2;
end

WAF = {};
WAF.Delay_resolution = Delay_resolution;
WAF.Doppler_resolution = Doppler_resolution;
WAF.Delay = Delay_WAF;
WAF.Doppler = Doppler_WAF;
WAF.WAF0 = WAF0;
WAF.WAF1 = WAF1;
WAF.Delta_Delay = Delta_Delay;
WAF.ACFbyACF = ACFbyACF;