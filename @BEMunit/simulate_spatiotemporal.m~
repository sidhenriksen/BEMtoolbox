function C = simulate_spatiotemporal(bem,generator,n_frames,duration,bootstrap_mode)
    % Computes the spatiotemporal response of the cell in response
    % to a trial consisting of the frames given by the stimulus
    % generator.

    dx = bem.x(2)-bem.x(1);
    dy = bem.y(2)-bem.y(1);

    if nargin < 3
        n_frames = 10;
    end

    if nargin < 4; 
        duration = 0.5;
    end

    if nargin < 5
        bootstrap_mode = false;
    end


    if bootstrap_mode
        % Then we need to figure out what we're doing here... If
        % we're in bootstrap mode, we need to 
        % 1) get the identifier for this BEMunit
        % 2) see if there is an existing match in the archive
        % 3a) if yes, load the match and use that data for the
        % Ls/Rs
        % 3b) if no, proceed as normal, but save the Ls and Rs, and
        % log the identifier + filename in the archive.
        
        id = bem.get_identifier();
        
        
        
    end

    n_frame_repeats = duration/(n_frames*bem.dt);           

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

        Ls = zeros(bem.n_subunits,n_frames*n_frame_repeats);
        Rs = zeros(bem.n_subunits,n_frames*n_frame_repeats);


        for k = 1:bem.n_subunits;
            % compute left and right eye responses
            L = bem.subunits(k).L(:)' * I_Ls;
            R = bem.subunits(k).R(:)' * I_Rs;

            L_rep = repmat(L,[n_frame_repeats,1]);
            R_rep = repmat(R,[n_frame_repeats,1]);

            Ls(k,:) = L_rep(:)*dx*dy;
            Rs(k,:) = R_rep(:)*dx*dy;
        end
    end
    
    
    
    V = zeros(size(Ls));
    
    for k = 1:bem.n_subunits;
        L = Ls(k,:);
        R = Rs(k,:);
        bem.t = linspace(0,(n_frames*n_frame_repeats-1)*bem.dt,n_frames*n_frame_repeats);
        rf = bem.subunits(k).rf_params;
        t_c = median(bem.t);
        switch bem.temporal_kernel
            case 'gaussian'
                L_tk=temporal_gaussian(bem.t,t_c,rf.left);
                R_tk=temporal_gaussian(bem.t,t_c,rf.right);

            case 'gamma-cosine'
                L_tk=gamma_cosine(bem.t,t_c,rf.left);
                R_tk=gamma_cosine(bem.t,t_c,rf.right);
                
            case 'none'
                error('Error: Assign a temporal kernel type (e.g. bem.temporal_kernel=''gaussian'')')
        end
        
        % Transform temporal kernel to frequency domain
        L_tk_fft = fft(L_tk);
        R_tk_fft = fft(R_tk);
        
        % Transform response to frequency domain
        L_fft = fft(L);
        R_fft = fft(R);
        
        % Multiply the ffts and move back to time domain; this is 
        % equivalent to convolving the two.
        V_L = ifft(L_tk_fft .* L_fft);
        V_R = ifft(R_tk_fft .* R_fft);
        
        % Finally, pass the summed response through the relevant
        % nonlinearity
        V(k,:) = bem.subunits(k).NL((V_L+V_R));        
    end
    C = sum(V);
end