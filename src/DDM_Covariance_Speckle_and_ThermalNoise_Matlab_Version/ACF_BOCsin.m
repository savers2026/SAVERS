function [CodePhase, ACF, sampling_rate, range_delay] = ACF_BOCsin(chip_rate, sampling_rate, length_in_chip, m, n)
    speedlight = 299792458.0;  % Spead of light, m/s
    chip_rate = 1.023e6*m;
    alpha = n/m;
    tc=1.0/alpha;
    code_phase = 2.0*(-length_in_chip:chip_rate/sampling_rate:length_in_chip);
    tmp = zeros(size(code_phase));
    tmp1 = zeros(size(code_phase));
    for k = 1-alpha:alpha-1
        for j = 1:length(code_phase)
            t=code_phase(j);
            tmp1(j) = 2.0*Tri(t*alpha-2.0*k,1.0)-Tri(t*alpha-2.0*k-1.0,1.0)-Tri(t*alpha-2.0*k+1.0,1.0);
        end
        tmp = tmp + tmp1*(alpha-abs(k));
    end

    CodePhase=code_phase/2/m;
    ACF=tmp/max(tmp);
    range_delay = CodePhase*speedlight/1.023e6; 
end

