 function C = simulate_spatial2(bem,generator,n_frames,bootstrap_mode,run_parallel)
    % Computes the spatial response of a BEMunit object to a sequence of
    % stimuli created the generator object.
    %
    % Usage: C = bem.simulate_spatial(generator,<n_frames>,<bootstrap_mode>,<seed_list>);    
    % generator: Stimulus generator object  
    % n_frames (optional): number of frames. Default is 1.    
    % bootstrap_mode (optional): Whether to save responses to disk to use for
    % resampling (typically used in conjunction with    
    % simulate_spatiotemporal). Default is 0 (off).   
    % seed_list (optional): A list of seeds used before creating the
    % stimuli. This can be useful if you need to recreate an exact stimulus
    % sequence.    

    if nargin < 3;
        n_frames = 1;
    end
    
    if nargin < 4
        bootstrap_mode = false;
    end
    
    if nargin < 5;
        run_parallel = false;
    end
    
    
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
            
            current_C = simulate_spatial(bem,generator,new_n_frames,bootstrap_mode);
                        
            C(start:stop) = current_C;
        end

    else
        

        if n_frames > 0
            S = zeros(bem.n_subunits,n_frames); % subunit responses
            L_idx = bem.subunits(1).L_active_idx;
            R_idx = bem.subunits(1).R_active_idx;
            I_Ls = zeros(length(L_idx),n_frames);
            I_Rs = zeros(length(R_idx),n_frames);
           

            if run_parallel
                parfor j = 1:n_frames;
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(L_idx);
                    I_Rs(:,j) = I_R(R_idx);
                end
            else
                
                for j = 1:n_frames;
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(L_idx);
                    I_Rs(:,j) = I_R(R_idx);
                end
            end
        end

        
        for k = 1:bem.n_subunits;
            
            if n_frames > 0
                % compute left and right eye responses
                L = bem.subunits(k).L(L_idx)' * I_Ls;
                R = bem.subunits(k).R(R_idx)' * I_Rs;
            elseif bootstrap_mode
                assert(isfield(bem.subunits(k),'V_L'),...
                    'Error: Bootstrap mode enabled with N=0, but no bootstrap distribution found.')
                L = bem.subunits(k).V_L;
                R = bem.subunits(k).V_R;
            else
                error('Error: Number of frames must be greater than 0');
            end

            if bootstrap_mode               
                bem.subunits(k).V_L = L;
                bem.subunits(k).V_R = R;
            end


            % add L and R, and pass through nonlinearity
            S(k,:) = bem.subunits(k).NL(L+R);                     
        end
        
        % Add together simple cell response to get complex response C
        C = sum(S);
    end    
    
    % create a dummy generator with corr=1 and dx=pref dx in order to 
    % scale responses
    dummy_generator = generator; 
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
    if isprop(bem,'outputNL');
        C = bem.outputNL(C)./norm_constant;
    else
        C = C./norm_constant;
    end
    
    
    
    if bootstrap_mode && n_frames > 0;
        % Save data to disk if bootstrap mode is enabled        
        bem.save_bootstrap(generator);                
    end
end