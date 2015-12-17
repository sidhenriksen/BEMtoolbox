function bem = load_bootstrap(bem,generator)
    % Method to load bootstrapped data into the BEMunit object, if
    % available.
    % Usage: bem = bem.load_bootstrap(generator);
    % generator: A stimulus generator object. 
    
    % Note: The generator object is saved in the 
    % bootstrap file in order to ensure that the saved values
    % were obtained from the same stimulus as the generator.
    % 

    bem_id = bem.get_identifier();
    stim_id = generator.get_identifier();

    % unique ID for bem + stim pairing
    idx_name = [num2str(bem_id),'.idx'];
    mat_name = [num2str(bem_id),'.mat'];

    idx_file = [bem.bootstrap_dir,'/',idx_name];
    mat_file = [bem.bootstrap_dir,'/',mat_name];
            
      
            
            
    % If there is already a file, we want to extend, not overwrite
    %fprintf('Saving file... ')
    if exist(mat_file,'file');
        % Load the index file
        fid = fopen(idx_file,'r');
        stim_index = textscan(fid,'%s'); 
        stim_index = stim_index{1};
        fclose(fid);
        
        index_match = strcmp(stim_index,num2str(stim_id));
        % If there is a match, load the .mat file
        if any(index_match);
            
            if ~bem.silent
                fprintf('File found in archive; loading samples.\n');
            end
            container = load(mat_file);
            big_bem = container.big_bem;
                        
            assert(isfield(big_bem,id2string(stim_id)),'Error: Generator ID not found in bootstrap file.');
            bstrap = big_bem.(id2string(stim_id));
            
            assert(bem.n_subunits == length(bstrap.bem.subunits),...
            'Error: Mismatching number of subunits. This is likely a bug.');
        
            for j = 1:bem.n_subunits
                bem_id = bem.get_identifier(j);
            
                match=0;
                for k = 1:length(bstrap.bem.subunits)
                    % might have to make this get_identifier(bstrap.bem,j)
                    bstrap_id = get_identifier(bstrap.bem,j);

                    if bem_id == bstrap_id;
                        match = j;
                        break
                    end

                end

                assert(match > 0, 'Error: No matching subunit found. This is likely a bug.');               

                bem.subunits(j).V_L = bstrap.bem.subunits(j).V_L;
                bem.subunits(j).V_R = bstrap.bem.subunits(j).V_R;            
            end

        end    
                        
    end

       
   
end