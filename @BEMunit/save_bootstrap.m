function save_bootstrap(bem,generator)
    % Saves the computed monocular responses (pre-temporal
    % filtering) to a file with a unique identifier for this model
    % + stimulus. Will extend the file if it already exists.,
    % Usage: bem.save_bootstrap(generator);
    % generator: a stimulus generator object

    bem_id = bem.get_identifier();
    stim_id = generator.get_identifier();

    % unique ID for bem + stim pairing
    idx_name = [num2str(bem_id),'.idx'];
    mat_name = [num2str(bem_id),'.mat'];

    idx_file = [bem.bootstrap_dir,'/',idx_name];
    mat_file = [bem.bootstrap_dir,'/',mat_name];


    % If there already exists a bem_index file then we want to append
    % it with the stim_id (but only if it doesn't already exist).
    % If it already exists then we shouldn't do anything.            
    if exist(idx_file,'file');
        fid = fopen(idx_file,'r');
        stim_index = textscan(fid,'%s'); 
        stim_index = stim_index{1};
        stim_match = strcmp(stim_index,num2str(stim_id));
        if ~any(stim_match);
            stim_index{length(stim_index)+1} = num2str(stim_id);
        end
    else
        stim_index = {num2str(stim_id)};
    end                        

    new_bem = bem;
    
    existing_bootstrap_file = bem.check_bootstrap(generator);
    % If there is an existing bootstrap file then we want to extend
    % it rather than override it.
    if existing_bootstrap_file
        if ~bem.silent
            fprintf('File found in archive; loading samples. ')
        end
        A = load(mat_file); big_bem = A.big_bem;
        bstrap_bem = big_bem.(id2string(stim_id)).bem;
        

        % merge the bstrap_bem and our current bem
        for k = 1:bem.n_subunits
            current_bem_id = get_identifier(bem,k);
            match = 0;
            for j = 1:length(bstrap_bem.subunits);
                bstrap_id = get_identifier(bstrap_bem,j);
                if bstrap_id == current_bem_id;
                    match = j;
                end

            end

            assert(match > 0,'No matching subunit found. This is likely a bug.');

            V_L = bem.subunits(k).V_L;
            V_R = bem.subunits(k).V_R;
            
            new_bem.subunits(k).V_L = [bstrap_bem.subunits(match).V_L,V_L];
            new_bem.subunits(k).V_R = [bstrap_bem.subunits(match).V_R,V_R];            
        end        
    else
        % If the model file already exists we don't want to overwrite it,
        % so we load it first and then add the details for the new
        % stimulus.
        if exist(mat_file,'file')
            load(mat_file);
        end                
    end

    big_bem.(id2string(stim_id)).generator = generator;
    big_bem.(id2string(stim_id)).bem = new_bem;

    if ~bem.silent
        fprintf('Saving file... ')
    end
    save(mat_file,'big_bem')

    fid=fopen(idx_file,'w');
    fprintf(fid,'%s\n',stim_index{:});
    fclose(fid);
    
    if ~bem.silent
        fprintf('Done.\n');
    end
            
end