function G = gamma_cosine(t,tc,rf)
    alpha = rf.alpha;
    omega = rf.omega;
    tau = rf.tau;
    t_phi = rf.t_phi;
    
    if isfield(rf,'t0');
        tc = rf.t0;
    end

        G = (1/(gamma(alpha)*tau^alpha) * (t-tc).^(alpha-1) ...
    .* exp(-(t-tc)/tau) .* cos(2*pi*omega*(t-tc) + t_phi));

    G = real(G);
    G(t<0) = 0;

end

