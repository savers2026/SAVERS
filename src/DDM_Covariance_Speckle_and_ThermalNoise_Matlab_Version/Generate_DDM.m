function DDM = Generate_DDM(WAF, PowerDD, DDM, thermal_noise_level)

%% ==============================================================================================
% To generate the DDM of the reflected signal by convolving the power in delay and Doppler domain and the WAF
% INPUT:
  % ACF: structure including the ACF parameters
  % PowerDD: reflected signal power in Delay and Doppler Domain
  % DDM: DDM parameters
  % thermal_noise_level: thermal noise level: kTb/Tc

% Output:
  % DDM: the delay Doppler map of the reflected signal
%=========================================================================================================

%% Power in delay and Doppler Domain and WAF
power_in_DD   = PowerDD.power_in_DD;
Doppler_Sigma = PowerDD.Doppler_Sigma;
delay_Sigma   = PowerDD.delay_Sigma;

WAF0 = WAF.WAF0.^2.0;
WAF1 = WAF.WAF1;

%% Normalized thermal noise level
noise_normalized = thermal_noise_level/sum(WAF0)/sum(WAF1);
% noise_normalized = thermal_noise_level;

%% 2-D convolution between power in DD and WAF
power_in_DD_plus_noise = (power_in_DD + noise_normalized);
% power_in_DD_plus_noise = (power_in_DD + noise_normalized)/sum(WAF0)/sum(WAF1);

DDM0 = zeros(size(power_in_DD_plus_noise));

for i = 1:length(Doppler_Sigma)
    DDM0(i,:) = conv(power_in_DD_plus_noise(i,:),WAF0,'same');
end
for i = 1:length(delay_Sigma)
    DDM0(:,i) = conv(DDM0(:,i),WAF1,'same');
end

DDM.meanDDM = DDM0(PowerDD.DDM_idx_Doppler, PowerDD.DDM_idx_delay);

