function out = quiescenceRateRampUp(is,quiescence_rate_max,oxygen_saturation,quiescence_threshold)

out = zeros(size(is));

max_rate_ind = is<=oxygen_saturation;
out(max_rate_ind) = quiescence_rate_max;
ramp_up_ind = ~max_rate_ind & is<quiescence_threshold;
out(ramp_up_ind) = quiescence_rate_max*(1 - (is(ramp_up_ind)-oxygen_saturation)/(quiescence_threshold-oxygen_saturation));