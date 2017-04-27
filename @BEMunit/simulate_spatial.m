 function C = simulate_spatial(bem,generator,n_frames,bootstrap_mode,run_parallel,seed)
    % C = bem.simulate_spatial(generator,<n_frames>,<bootstrap_mode>,<run_parallel>,<seed_list>);    
    %
    % Computes the spatial response of a BEMunit object to a sequence of
    % stimuli created the generator object.
    % generator: Stimulus generator object  
    % n_frames (optional): number of frames. Default is 1.    
    % bootstrap_mode (optional): Whether to save responses to disk to use for
    % resampling (typically used in conjunction with    
    % simulate_spatiotemporal). Default is 0 (off).   
    % seed_list (optional): A list of seeds used before creating the
    % stimuli. This can be useful if you need to recreate an exact stimulus
    % sequence.
    
    if isnumeric(generator);
        S = zeros(size(generator,1),1);
        for k = 1:length(bem.subunits);
            
            B = [bem.subunits(k).L(:);bem.subunits(k).R(:)];
            
            S = S+bem.subunits(k).NL(generator*B);
           
        end        
        
        C = S;
        return;
    end

    assert((bem.Nx == generator.Nx) && (bem.Ny == generator.Ny),'Error: BEMunit and generator object dimensions must be equal.');
    
    if nargin < 3;
        n_frames = 1;
    end
    
    if nargin < 4
        bootstrap_mode = false;
    end
    
    if nargin < 5;
        run_parallel = false;
    end
    
    if nargin < 6
        seed=randi(1e9,1);        
    end    
    rng(seed)    
    seed_list = randi(1e9,n_frames,1);
    
    total_bytes = bem.Nx * bem.Ny *n_frames * 8;
    if total_bytes > bem.memory_threshold;
        % In order not flood memory, we split this into smaller pieces and
        % recursively call simulate_spatial with n_frames reduced.
        % We then iteratively place those in the response vector C
        % before exiting.
        
        C = zeros(1,n_frames);
        
        if ~bem.silent
            warning(['Dimensionality of stimulus too high; using serialised implementation to conserve memory (this is slower). ', ...
             'You can manually increase the limit by changing bem.memory_threshold property.'])
        end
        
        repeats = ceil(2*total_bytes/bem.memory_threshold);
        new_n_frames = round(n_frames/(2*total_bytes)*bem.memory_threshold);
        for k = 1:repeats;
            
            if k == repeats;
                new_n_frames = n_frames-(k-1)*new_n_frames;
            end
            
            start = (k-1)*new_n_frames + 1;
            stop = start+new_n_frames-1;
            
            
            current_C = simulate_spatial(bem,generator,new_n_frames,bootstrap_mode,-1);
                        
            C(start:stop) = current_C;

        end
        % The bootstrap data here is saved to disk in batches and is
        % not available in the workspace, so we have to load it from
        % disk.
        if bootstrap_mode
            bem = bem.load_bootstrap(generator); 
        end
    else
        

        if n_frames > 0
            S = zeros(bem.n_subunits,n_frames); % subunit responses
            Ls = zeros(bem.n_subunits,n_frames);
            Rs = zeros(bem.n_subunits,n_frames);
            I_Ls = zeros(bem.Nx*bem.Ny,n_frames);
            I_Rs = zeros(bem.Nx*bem.Ny,n_frames);

            if run_parallel==1
                parfor j = 1:n_frames;
                    
                    rng(seed_list(j));
                    
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(:);
                    I_Rs(:,j) = I_R(:);
                end
            else
                
                for j = 1:n_frames;
                    rng(seed_list(j));
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(:);
                    I_Rs(:,j) = I_R(:);
                end
            end
        end

        
        for k = 1:bem.n_subunits;
            
            if n_frames > 0
                % compute left and right eye responses
                L = bem.subunits(k).L(:)' * I_Ls;
                R = bem.subunits(k).R(:)' * I_Rs;
            elseif bootstrap_mode
                assert(isfield(bem.subunits(k),'V_L'),...
                    'Error: Bootstrap mode enabled with N=0, but no bootstrap distribution found.')
                L = bem.subunits(k).V_L;
                R = bem.subunits(k).V_R;
            else
                error('Error: Number of frames must be greater than 0');
            end

            if bootstrap_mode               
                bem.subunits(k).V_L = L';
                bem.subunits(k).V_R = R';
            end
          
            % add L and R, and pass through nonlinearity
            Ls(k,:) = L;
            Rs(k,:) = R;
            
        end
        
        if bem.norm_params(1)
            norm_Ls = zeros(size(Ls));
            norm_Rs = zeros(size(Rs));
            for k = 1:bem.n_subunits;            
                % Find the appropriate indices for this cell (to make a
                % monocular complex cell)
                idx = find([bem.quad_pairs{:}] == k);
                current_idx = floor((idx+1)/2);
                current_pair = bem.quad_pairs{current_idx};

                
                L_complex = sum(Ls(current_pair,:).^2);
                R_complex = sum(Rs(current_pair,:).^2);                

                L_norm = 1./(1+bem.norm_params(1)*L_complex.^bem.norm_params(2));
                R_norm = 1./(1+bem.norm_params(1)*R_complex.^bem.norm_params(2));

                norm_Ls(k,:) = Ls(k,:).*L_norm;
                norm_Rs(k,:) = Rs(k,:).*R_norm;
                
            end
            Ls = norm_Ls;
            Rs = norm_Rs;
        end
    
    
        for k = 1:bem.n_subunits;
            S(k,:) = bem.subunits(k).NL(Ls(k,:)+Rs(k,:))*bem.subunits(k).weight;
        end

    
        % Add together simple cell response to get complex response C
        if bem.n_subunits > 1
            C = sum(S);
        else
            C = S;
        end
    end    
    
    % create a dummy generator with corr=1 and dx=pref dx in order to 
    % scale responses
    dummy_generator = generator.copy(); 
    dummy_generator.correlation = 1;
    dummy_generator.dx = round(bem.dx/bem.deg_per_pixel);
    % This will see if there is a pre-computed normalisation term.
    % If there isn't, then we will just set this to 1.
    idstring = id2string(dummy_generator.get_identifier());    
    
    norm_constant = 1;
    if isprop(bem,'norm_constants');
        if isfield(bem.norm_constants,idstring);
            norm_constant = bem.norm_constants.(idstring);              
        end
    end
    
    % Pass through output nonlinearity before returning
    % Have to do this because of some backwards compatibility issues
    if isprop(bem,'outputNL')
        C = bem.outputNL(C)./norm_constant;
    else
        C = C./norm_constant;
    end
    
    
    if bootstrap_mode && n_frames > 0
        % Save data to disk if bootstrap mode is enabled and we didn't
        % serialise it
        
        % So this will be true for either recursive calls (run_parallel ==-1)
        % OR for method calls that never made any recursve calls (bytes <= thresh).
        if run_parallel == -1 || (total_bytes <= bem.memory_threshold)
            bem.save_bootstrap(generator);
        end
    end
end