function clean_bootstrap(bem)

    bem_id = bem.get_identifier();

    % unique ID for bem + stim pairing
    idx_name = [num2str(bem_id),'.idx'];
    mat_name = [num2str(bem_id),'.mat'];

    idx_file = [bem.bootstrap_dir,idx_name];
    mat_file = [bem.bootstrap_dir,mat_name];
            
    changed = 0;
            
    % If there is already a file, we want to extend, not overwrite
    %fprintf('Saving file... ')
    if exist(mat_file,'file');
        
        
        fid = fopen(idx_file,'r');
        stim_index = textscan(fid,'%s'); 
        stim_index = stim_index{1};
        fclose(fid);
        
        container = load(mat_file);
        big_bem = container.big_bem;
        
        bem_fields = fields(big_bem);
        
        for j = 1:length(bem_fields);
            
            current = big_bem.(bem_fields{j});
            
            if isempty(current.bem);
                big_bem=rmfield(big_bem,bem_fields{j});
                
                
                
                if ~bem.silent
                    fprintf('Found incomplete ID; deleting.\n');
                end
                
                idx = ~strcmp(stim_index,num2str(string2id(bem_fields{j})));
                stim_index = stim_index(idx);
                changed = 1;
            end
            
        end
        
    end
    
    
    if changed
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

end