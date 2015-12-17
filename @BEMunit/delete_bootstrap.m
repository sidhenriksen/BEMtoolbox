function bem = delete_bootstrap(bem)
    % Method to delete bootstrap data for this file.
    % This can be useful if the data is corrupt. Use with caution,
    % once the data is deleted it cannot be recovered.
    % Usage: bem = bem.delete_bootstrap()
        
    
    bem_id = bem.get_identifier();
    
    % unique ID for bem + stim pairing
    idx_name = [num2str(bem_id),'.idx'];
    mat_name = [num2str(bem_id),'.mat'];

    idx_file = [bem.bootstrap_dir,'/',idx_name];
    mat_file = [bem.bootstrap_dir,'/',mat_name];
            
    % If there is already a file, we want to extend, not overwrite
    %fprintf('Saving file... ')
    if exist(mat_file,'file');
        delete(mat_file);
        delete(idx_file);
        fprintf ('Data deleted successfully.\n');
    else
        fprintf('No data found. Exiting...\n');
       
    end
   
end