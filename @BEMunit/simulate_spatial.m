 function C = simulate_spatial(bem,generator,n_frames)

    if nargin < 3;
        n_frames = 1;
    end


    S = zeros(bem.n_subunits,n_frames); % subunit responses

    total_bytes = bem.Nx * bem.Ny *n_frames * 8;
    if total_bytes > bem.memory_threshold;
        warning(['Dimensionality of stimulus too high; using serialised implementation to conserve memory (this is slower). ', ...
            'You can manually increase the limit by changing bem.memory_threshold property.'])
        error('Serialisation not yet implemented. Reduce dimensionality or increase memory threshold');

    else

        I_Ls = zeros(bem.Nx*bem.Ny,n_frames);
        I_Rs = zeros(bem.Nx*bem.Ny,n_frames);

        for j = 1:n_frames;
            [I_L,I_R] = generator.generate();
            I_Ls(:,j) = I_L(:);
            I_Rs(:,j) = I_R(:);
        end    

        for k = 1:bem.n_subunits;
            % compute left and right eye responses
            L = bem.subunits(k).L(:)' * I_Ls;
            R = bem.subunits(k).R(:)' * I_Rs;

            % add L and R, and pass through nonlinearity
            S(k,:) = bem.subunits(k).NL(L+R);                     
        end
    end

    % Add together subunit responses to get final response
    % (complex cell response in the case of an energy model unit)
    C = sum(S);

end