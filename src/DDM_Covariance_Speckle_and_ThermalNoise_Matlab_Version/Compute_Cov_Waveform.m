function Cov_Z_WF = Compute_Cov_Waveform(WAF, PowerDD, DDM, thermal_noise_level, Tc, Tinc, freq, freq_prime)

%% ==============================================================================================
% To compute the covarance of one waveform
% The covariance of the DDM is:
  % Cov_DDM(tau_0, tau_1, f_0, f_1) = <DDM(tau_0, f_0)*DDM(tau_1, f_1)>
% Given the frequencies f_0 and f_1, this function is to compute:
  % Cov_DDM(tau_0, tau_1)
% INPUT:
  % WAF: structure including the ACF parameters
  % PowerDD: reflected signal power in Delay and Doppler Domain
  % DDM: DDM parameters
  % thermal_noise_level: thermal noise level: kTb/Tc
  % Tc: coherent integration time
  % Tinc: ioncoherent average time
  % freq: frequency 1
  % freq_prime: frequency 2
% Output:
  % Cov_Z_WF: covariance of the waveform, an array with N_delay_bin * N_delay_bin elements
%=========================================================================================================

%% number of delay bins
    delay_resolution = DDM.Delay_resolution;
    N_delay_bin = length(DDM.DDM_delay);
    Cov_Z_WF = zeros(N_delay_bin, N_delay_bin);
    Delta_Delay = (0:N_delay_bin-1)*delay_resolution;

%% Compute the waveform covariance by diagonals
    for d_delay = Delta_Delay
        delta_idx = round(d_delay/delay_resolution);
        if d_delay < 400.0
            Cov_Z = Compute_Cov_diagonal(WAF, PowerDD, thermal_noise_level, Tc, Tinc, d_delay, freq, freq_prime);
        else
            Cov_Z = zeros(1, 301);
        end
        len_diag = N_delay_bin - delta_idx;
        if delta_idx == 0
            Cov_Z_WF = Cov_Z_WF + diag(Cov_Z(1:len_diag),delta_idx);
        else
            Cov_Z_WF = Cov_Z_WF + diag(Cov_Z(1:len_diag),delta_idx) + diag(Cov_Z(1:len_diag),-delta_idx);
        end
    end