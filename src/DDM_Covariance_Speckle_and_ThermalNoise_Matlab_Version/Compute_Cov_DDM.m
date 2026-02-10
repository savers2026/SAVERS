function Cov_Z_DDM = Compute_Cov_DDM(WAF, PowerDD, DDM, thermal_noise_level, Tc, Tinc)

%% ==============================================================================================
% To compute the covarance of the DDM
% The covariance of the DDM is:
  % Cov_DDM(tau_0, tau_1, f_0, f_1) = <DDM(tau_0, f_0)*DDM(tau_1, f_1)>
% INPUT:
  % WAF: structure including the ACF parameters
  % PowerDD: reflected signal power in Delay and Doppler Domain
  % DDM: DDM parameters
  % thermal_noise_level: thermal noise level: kTb/Tc
  % Tc: coherent integration time
  % Tinc: ioncoherent average time
% Output:
  % Cov_Z_DDM: covariance of the waveform, an 4-D array with N_delay_bin*N_Doppler_bin*N_delay_bin*N_Doppler_bin elements
%=========================================================================================================

%% Size of the DDM covariance and initializations
    N_delay_bin = length(DDM.DDM_delay);
    N_Doppler_bin = length(DDM.DDM_Doppler);
    DDM_Doppler = DDM.DDM_Doppler;
    Cov_Z_DDM = zeros(N_delay_bin, N_Doppler_bin, N_delay_bin, N_Doppler_bin);

%% Compute the DDM covariance by waveform at each Doppler bin
    for i = 1:N_Doppler_bin
        for j = i:N_Doppler_bin
            Cov_Z_DDM(:,i,:,j) = Compute_Cov_Waveform(WAF, PowerDD, DDM, thermal_noise_level, Tc, Tinc, DDM_Doppler(i), DDM_Doppler(j));
            Cov_Z_DDM(:,j,:,i) = Cov_Z_DDM(:,i,:,j);
        end
    end
end
