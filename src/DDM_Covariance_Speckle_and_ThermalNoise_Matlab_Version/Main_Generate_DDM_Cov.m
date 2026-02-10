
%% --------------------------------------------------------------------------------------
% -------------------------------Constants---------------------------------------------
%--------------------------------------------------------------------------------------

speedlight        = 299792458.0;  % Spead of light, m/s
%%%%LINKED - Double Check this
chip_rate         = 1/(geoSYSp.chipTime*1e-6);%%%1.023e6; %E1, L1
%%%%to be LINKED ??? 
sampling_rate_acf = 1.023e8; %E1, L1
%%%%to be LINKED ??? 
length_in_chip    = 2.0;

BOC_m             = 1;
BOC_n             = 1;


%% --------------------------------------------------------------------------------------
% -------------------------------Input Variables---------------------------------------
% ---------------Either system parameters or the outputs from SAVERS-------------------
%--------------------------------------------------------------------------------------

% Instrument parameters: SYSTEMATIC PARAMETERS


%%%%to be LINKED ??? 
rcv_bw = 4.0e6;   % Receiver bandwidth: 4.0 MHz
%%%% VARIABLE LINKED
Tc     = geoSYSp.cohIntTime*1e-3;%%0.001;   % Coherent Integration Time: 0.001 second
Tinc   = geoSYSp.incoIntTime*1e-3;%%%1.0;   % Incoherent averaging time: 1.0 second

% DDM parameters: SYSTEMATIC PARAMETERS
%%%%to be LINKED ??? 
% DDM delay sampling rate 16.367 MHz (TDS-1 configuration) 
%Delay_resolution = speedlight/16.367e6;  
%HydroGNSS is 53 MHz => 18.87 ns => 5.661 metres. (Hydro configuration) 
% Delay_resolution = speedlight/53e6;
% Delay_resolution =  0.5*speedlight*(DDMparam.fdelay(1,1)- DDMparam.fdelay(1,2))*1e-6;
Delay_resolution = 73;
% DDM delay bins index, from 1 to 128, the 31st delay bin corresponds to the specular point
% delay_bins       = (1:128) - 31;   

%%%QUI VALE IN GENERALE?
delay_bins      = (1:128) - 65;   
DDM_delay        = delay_bins*Delay_resolution;

% Rplot
% Fplot

% geoSYSp.nx_Rdelay
% geoSYSp.ny_Dfreq

%DDM Doppler resolution, 200.0 Hz
%Doppler_resolution = 500.0;
Doppler_resolution = geoSYSp.fpas;
%DDM Doppler bins index, from 1 to 21, the 11st Doppler bin correspoinds to
%the specular point
Doppler_bins       = (1:20) - 11;  
%Doppler_bins      = (1:20) - 0;  
DDM_Doppler        = Doppler_bins*Doppler_resolution;

% To put DDM parameters into a structure
DDM = {};
DDM.Delay_resolution   = Delay_resolution;
DDM.delay_bins         = delay_bins;
DDM.DDM_delay          = DDM_delay;
DDM.Doppler_resolution = Doppler_resolution;
DDM.Doppler_bins       = Doppler_bins;
DDM.DDM_Doppler        = DDM_Doppler;

% Delay, Doppler and Reflected signal power of each grid: OUTPUTS LINKED

%%%%LINKED - Double Check this
Relative_Doppler_grid = DDMparam.fdoppler;
delay_grid            = speedlight*DDMparam.fdelay*1e-6; %in meters

%%%%LINKED - Double Check this
% sigPower_noNoise_LR_lin
%ref_power_grid                = sigmaLR_cohe + sigmaLR_inco;
% ref_power_grid               = sigPower_noNoise_LR_lin;
ref_power_grid                 = DDMparam.PonSurf;
 
Refl_Surf                      = {};
Refl_Surf.delay_surface_grid   = delay_grid;
Refl_Surf.Doppler_surface_grid = Relative_Doppler_grid;
% Refl_Surf.power_surface_grid   = ref_power_grid*Tc;
Refl_Surf.power_surface_grid   = ref_power_grid*Tc*2;

% Thermal noise level: OUTPUTS FROM SAVERS
%5.9354e-18; % Typical thermal noise level, Watts, k*Tb/Tc
%%%%%% To BE LINKED
thermal_noise_level = 0;
%%%setting thermal_noise_level = 0 we still get a Covariance map, which models
%%%the 'speckle' in the variable Refl_Surf


%% -------------------------------------------------------------------------------------
% Generate the ACF (Ideal and Bandlimited) for GPS CA code and Galileo E1B/C signals
%--------------------------------------------------------------------------------------
[code_phase, ACF_CA, sampling_rate_acf, ACF_CA_range_delay]   = ACF_BPSK(chip_rate, sampling_rate_acf, length_in_chip);                 %%L1
%%to do E5-L5 just change chip_rate. For E5 = 1/0.1e-6 = 10e6;
[code_phase, ACF_BOC, sampling_rate_acf, ACF_BOC_range_delay] = ACF_BOCsin(chip_rate, sampling_rate_acf, length_in_chip, BOC_m, BOC_n); %%E1


% Bandlimited ACFs
%[b,a] = butter(6, rcv_bw/(sampling_rate_acf/2));
%ACF_CA_BW = filter(b, a, ACF_CA);
%ACF_BOC_BW = filter(b, a, ACF_BOC);
% To put ACF parameters into a structure
ACF = {};
ACF.code_phase        = code_phase;
ACF.ACF               = ACF_CA;
ACF.sampling_rate_acf = sampling_rate_acf;
ACF.ACF_range_delay   = ACF_CA_range_delay;


%% -------------------------------------------------------------------------------------
% Reflected power in Delay and Doppler Domain, WITHOUT convolution with the WAF
%--------------------------------------------------------------------------------------
%%% Refl_Surf is defined over an xy grid (on the surface) and it assumes
%%% the radar equation has been already applied
[PowerDD, WAF] = MapinDDandWAF(ACF, Refl_Surf, DDM, Tc);
%%%PowerDD is defined over the D and d map (grid). So that the grid has a
%%%different size wrt Refl_Surf

%PowerDD.power_in_DD = PowerDD.power_in_DD./100;

%% -------------------------------------------------------------------------------------
% Power DDM computation
%--------------------------------------------------------------------------------------
DDM = Generate_DDM(WAF, PowerDD, DDM, thermal_noise_level);

%DDM.meanDDM = ones(21,128);

%% -------------------------------------------------------------------------------------
% Generate the covariance of the DDM -- signal + thermal noise
%--------------------------------------------------------------------------------------
Cov_Z_DDM = Compute_Cov_DDM(WAF, PowerDD, DDM, thermal_noise_level, Tc, Tinc);
%Cov_Z_WF = Compute_Cov_Waveform(WAF, PowerDD, DDM, thermal_noise_level, Tc, Tinc, 0.0, 0.0);
Cov_Z_DDM_flatten = reshape(Cov_Z_DDM,[length(DDM.delay_bins)*length(DDM.DDM_Doppler),length(DDM.delay_bins)*length(DDM.DDM_Doppler)]);


%% Test added by DC

%figure, imagesc(Cov_Z_DDM(:,:,1,1)), colorbar
%figure, imagesc(Cov_Z_DDM_flatten), colorbar

Fplot = 1:1:geoSYSp.ny_Dfreq;
Rplot = 1:1:geoSYSp.nx_Rdelay;

% DDM.meanDDM = DDM.meanDDM/(500*500);

%sometime there is a very small Im - the plot fails
% DDM.meanDDM = real(DDM.meanDDM);

%sometime there is a negative value - the dB plot fails
DDM.meanDDM(DDM.meanDDM < 0) = nan;

figure, imagesc(Rplot, Fplot, 10*log10(DDM.meanDDM)), colorbar, axis xy
caxis([max(10*log10(DDM.meanDDM(:))) - 20, max(10*log10(DDM.meanDDM(:)))])
colormap(jet)
title('mean DDM')
%%% check the different noise floor in DDM.meanDDM

figure, imagesc(Refl_Surf.power_surface_grid), colorbar, axis xy
title('power on surface')
% figure, imagesc(PowerDD.power_in_DD), colorbar, axis xy


%%%next step: output with a noisy DDM

%%% Noisy DDM - added in a second phase

Cov_DDM = (Cov_Z_DDM_flatten + Cov_Z_DDM_flatten') / 2;
[V,D]   = eig(Cov_DDM);
D1      = abs(D);

Cov_DDM           = V*D1*V';
DDM_mean          = DDM.meanDDM';
mean_ddm          = DDM_mean(:);
noisy_DDM_flatten = mvnrnd(mean_ddm,Cov_DDM,1);
noisy_DDM         = reshape(noisy_DDM_flatten, size(DDM.meanDDM'));

%figure, imagesc(Doppler_bins, delay_bins, noisy_DDM.'), colorbar
figure, imagesc(Rplot, Fplot, noisy_DDM.'), colorbar
%caxis([max(10*log10(DDM.meanDDM(:))) - 20, max(10*log10(DDM.meanDDM(:)))])
colormap(jet)
title('noisy DDM')

