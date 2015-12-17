function C = simulate_spatiotemporal(bem,generator,n_frames,duration,bootstrap_mode)
    % Computes the spatiotemporal response of the cell in response
    % to a trial consisting of the frames given by the stimulus
    % generator.
    % Usage: C = bem.simulate_spatiotemporal(generator,<n_frames>,<duration>,<bootstrap_mode>,<seed_list>)
    % generator: stimulus generator object (e.g. RDS, pairedRDS, etc.)
    % n_frames (optional): number of frames. Default is 10.
    % duration (optional): duration of stimulus. Default is 500ms.
    % bootstrap_mode (optional): whether to compute responses using
    % resampling (requires to run simulate_spatial in bootstrap mode first)
    
    

    dx = bem.x(2)-bem.x(1);
    dy = bem.y(2)-bem.y(1);

    if nargin < 3
        n_frames = 10;
    elseif isempty(n_frames);
        n_frames = 10;
    end

    if nargin < 4; 
        duration = 0.5;
    elseif isempty(duration);
        duration = 0.5;
    end

    if nargin < 5
        bootstrap_mode = false;
    end
    
    n_frame_repeats = ceil(duration/(n_frames*bem.dt));

    total_bytes = bem.Nx * bem.Ny *n_frames * 8;

    if total_bytes > bem.memory_threshold;
        warning(['Dimensionality of stimulus too high; using serialised implementation to conserve memory (this is slower). ', ...
            'You can manually increase the limit by changing bem.memory_threshold property.'])
        error('Serialisation not yet implemented. Reduce dimensionality or increase memory threshold');

    else

        
        Ls = zeros(bem.n_subunits,n_frames*n_frame_repeats);
        Rs = zeros(bem.n_subunits,n_frames*n_frame_repeats);

        % If we're not doing bootstrap, we're just gonna run this as normal
        if ~bootstrap_mode
            
            % First generate some n_frames images
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

                L_rep = repmat(L,[n_frame_repeats,1]);
                R_rep = repmat(R,[n_frame_repeats,1]);

                Ls(k,:) = L_rep(:)*dx*dy;
                Rs(k,:) = R_rep(:)*dx*dy;
            end
        else
            assert(isfield(bem.subunits,'V_L'),'Error: Bootstrapped data not loaded.');                        
            seed_indices = randi(length(bem.subunits(1).V_L),1,n_frames);
            
            % If we're in bootstrap mode and we've specified a seed
            % sequence then we find the indices for those and use them as
            % the responses

            for k = 1:bem.n_subunits;
                % The seed indices here are either given as arguments, or
                % if no argument is given, are random picks from
                % V_L/V_R.
                L = bem.subunits(k).V_L(seed_indices);
                R = bem.subunits(k).V_R(seed_indices);
                L_rep = repmat(L,[n_frame_repeats,1]);
                R_rep = repmat(R,[n_frame_repeats,1]);

                Ls(k,:) = L_rep(:)*dx*dy;
                Rs(k,:) = R_rep(:)*dx*dy;
            end
        end
            
    end
    
    
    
    V = zeros(size(Ls));
    
    if strcmp(bem.temporal_kernel,'none')
        warning('No temporal kernel assigned; defaulting to Gaussian.');
        bem.temporal_kernel = 'gaussian';
    end
    
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
        end
        
        % Transform temporal kernel to frequency domain
        L_tk_fft = fft(L_tk);
        R_tk_fft = fft(R_tk);
        
        % Transform response to frequency domain
        L_fft = fft(L);
        R_fft = fft(R);
        
        % Multiply the ffts and move back to time domain; this is 
        % equivalent to convolving the two.
        V_L = fftshift(ifft(L_tk_fft .* L_fft));
        V_R = fftshift(ifft(R_tk_fft .* R_fft));
        
        % Finally, pass the summed response through the relevant
        % nonlinearity
        V(k,:) = bem.subunits(k).NL((V_L+V_R));        
    end
    
    C = sum(V);
    
    % This will check if there is a pre-computed normalisation term.
    % If there isn't, then we will just set this to 1.
    
    
    if isprop(generator,'correlation') && isprop(generator,'dx')
        dummy_generator = generator; 
        dummy_generator.correlation = 1;
        dummy_generator.dx = round(bem.dx/bem.deg_per_pixel);
        
        idstring = id2string(dummy_generator.get_identifier());
    
        if isfield(bem.norm_constants,idstring);
            norm_constant = bem.norm_constants.(idstring);
        else
            norm_constant = 1;
        end
    else
        norm_constant = 1;
    end
    
    C = bem.outputNL(C)/norm_constant;
    C = C(1:round(duration/bem.dt));
end