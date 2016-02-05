function G = temporal_gaussian(t,rf)

    G = exp(-(t-rf.t0).^2 ./ (2*rf.tau^2));
end