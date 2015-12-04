function G = gamma_cosine(t,tc,rf)
    alpha = rf.alpha;
    omega = rf.omega;
    tau = rf.tau;
    t_phi = rf.t_phi;

        G = (1/(gamma(alpha)*tau^alpha) * (t-tc).^(alpha-1) ...
    .* exp(-(t-tc)/tau) .* cos(2*pi*omega*(t-tc) + t_phi));

    G = real(G); 
    G(t<tc) = 0; 

end

