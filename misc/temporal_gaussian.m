function G = temporal_gaussian(t,tc,rf)
    tau = rf.tau;

    G = exp(-(t-tc).^2 ./ (2*tau^2));
end