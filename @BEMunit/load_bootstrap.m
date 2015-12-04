function bem = load_bootstrap(bem)
    % Method to load bootstrapped data into the BEMunit object, if
    % available.
    % Usage: bem = bem.load_bootstrap();
    % 

    files = what(bem.bootstrap_dir);
    mat_files = files.mat;
    
    id = bem.get_identifier();
    
    file_match = strcmp(mat_files,[id,'.mat']);
    
    if any(file_match);
        fprintf('File found in archive; loading samples.\n');
        
        bootstrap_file = load(mat_files(file_match));
    else
        fprintf('No match found. Exiting.\n');
    end
end