function G = gamma_cosine(t,rf)
    alpha = rf.alpha;
    omega = rf.omega;
    tau = rf.tau;
    t_phi = rf.t_phi;
    
    if isfield(rf,'t0');
        rf.t0 = rf.t0;
    end

        G = (1/(gamma(alpha)*tau^alpha) * (t-rf.t0).^(alpha-1) ...
    .* exp(-(t-rf.t0)/tau) .* cos(2*pi*omega*(t-rf.t0) + t_phi));

    G = real(G);
    G(t<0) = 0;

end

