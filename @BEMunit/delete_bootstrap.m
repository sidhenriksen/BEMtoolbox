function delete_bootstrap(bem,generator)
    % Method to delete bootstrap data for this file.
    % This can be useful if the data is corrupt. Use with caution,
    % once the data is deleted it cannot be recovered.
    % Usage: bem = bem.delete_bootstrap(<generator>)
    % generator: optional argument; will only delete for this generator   
    
    bem_id = bem.get_identifier();
    
    if nargin > 1
        stim_id = generator.get_identifier();
        
        csvfiles = {[bem.bootstrap_dir,num2str(bem_id),'_',num2str(stim_id),'.csv']};
    else        
        csvfiles = dir([bem.bootstrap_dir,num2str(bem_id),'_*.csv']);
    end
    
        
    if ~isempty(csvfiles)

        if length(csvfiles)==1;
            if iscell(csvfiles)
                delete(csvfiles{1});
            else
                delete([bem.bootstrap_dir,csvfiles.name])
            end
            fprintf('1 file deleted successfully.\n');
        else
            for j = 1:length(csvfiles);            
                delete([bem.bootstrap_dir,csvfiles(j).name]);
            end        
            fprintf('%i files deleted successfully.\n',length(csvfiles));
        end
        
    else
        fprintf('No data found. Exiting...\n');       
    end
   
end