function Cov_Z = Compute_Cov_diagonal(WAF, PowerDD, thermal_noise_level, Tc, Tinc, d_delay, freq, freq_prime)

%% ==============================================================================================
% To generate the covariance of the DDM along the diagonal
% The covariance of the DDM is:
  % Cov_DDM(tau_0, tau_1, f_0, f_1) = <DDM(tau_0, f_0)*DDM(tau_1, f_1)>
% Given the frequencies f_0 and f_1, this function is to compute the elements along each diagonal, i.e.
  % Cov_DDM(tau_0 - tau_1, f_0, f_1)
% INPUT:
  % WAF: structure including the WAF parameters
  % PowerDD: reflected signal power in Delay and Doppler Domain
  % thermal_noise_level: thermal noise level: kTb/Tc
  % Tinc: incoherent average time
  % Tc: coherent integration time
  % d_delay: delay difference between the waveform lags
  % freq, freq_prime: frequencies for the computation of the DDM covariance
% Output:
  % WAF: Woodward's ambiguity function adapted to the DDM resolution
  % PowerDD: reflected signal power in Delay and Doppler Domain
%=========================================================================================================

%% y: complex DDM; z: one-shot power waveform, Z: averaged power waveform
%% To compute ACF(tau)*ACF(tau + delta_tau)
chiplength = 299792458.0/1.023e6;
delay_resolution_itp = PowerDD.delay_resolution_itp;
power_in_DD = PowerDD.power_in_DD;
Doppler_Sigma = PowerDD.Doppler_Sigma;
delay_resolution = WAF.Delay_resolution;

n_sample = ceil(1.5*chiplength/delay_resolution_itp);
Delay = (-n_sample:n_sample)*delay_resolution_itp;
ACFbyACF = interp1(WAF.Delay, WAF.ACFbyACF(round(d_delay/delay_resolution)+1,:), Delay, 'PCHIP', 0);


%% Sigma_cov_ACFbyACF is the convolution between ACFbyACF and power in Delay-Doppler Domain [Sigma(delay, f)] --signal component
Sigma_cov_ACFbyACF = zeros(size(power_in_DD));
for i = 1:length(Doppler_Sigma)
    Sigma_cov_ACFbyACF(i,:) = conv(power_in_DD(i,:),ACFbyACF,'same');
end
Sigma_cov_ACFbyACF = Sigma_cov_ACFbyACF(:, PowerDD.DDM_idx_delay);

%% Sigma_cov_ACFbyACF is the convolution between ACFbyACF and power in Delay-Doppler Domain [Sigma(delay, f)] --thermal noise component
N0_by_ACF = interp1(WAF.Delay, WAF.WAF0, d_delay, 'linear', 0.0)*thermal_noise_level/sum(WAF.WAF1);

%% Thermal noise and signal
Sigma_cov_ACFbyACF = Sigma_cov_ACFbyACF + N0_by_ACF;

%%  SincbySinc = Sinc(f)*Sinc(f + delta_f)
SincbySinc = sinc(Tc*(Doppler_Sigma + freq)).*sinc(Tc*(Doppler_Sigma + freq_prime));
SincbySinc1 = repmat(SincbySinc, [length(PowerDD.DDM_idx_delay), 1]);

%% covariance of the complex DDM in frequency domain: C_y(xi, d_delay, freq, freq_prime)
Cov_y_freq_domain = Sigma_cov_ACFbyACF.*(SincbySinc1');

%% covariance of the complex DDM in time domain: C_y(d_t, d_delay, d_freq)
Cov_y_time_domain = abs(fftshift(fft(Cov_y_freq_domain, [], 1), 1));
freq_range = Doppler_Sigma(end) - Doppler_Sigma(1);

%% Cov_z: covariance of one-shot power DDM
Cov_z = Cov_y_time_domain.^2.0;
      
%% Cov_Z: covariance of the averaged power DDM
N_inc = Tinc/Tc;
Cov_Z = sum(Cov_z(1:(Tc*freq_range):end,:),1)/N_inc;
