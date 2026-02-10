function [code_phase, Lambda, sampling_rate, range_delay] = ACF_BPSK(chip_rate, sampling_rate, length_in_chip)
    speedlight = 299792458.0;  % Spead of light, m/s
    code_phase = -length_in_chip:chip_rate/sampling_rate:length_in_chip;
    Lambda = 1.0 -abs(code_phase);
    Lambda(abs(code_phase)>1.0) = 0.0;
    range_delay = code_phase*speedlight/chip_rate;
end
