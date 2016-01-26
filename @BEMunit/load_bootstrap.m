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
    
    csvfile = [bem.bootstrap_dir,num2str(bem_id),'_',num2str(stim_id),'.csv'];      
    
    if exist(csvfile,'file');
        if ~bem.silent
            fprintf('Bootstrap data found; loading.\n');
        end
        
        data = csvread(csvfile);
    
        for j = 1:length(bem.subunits);
            k = (j-1)*2 + 1;
            
            bem.subunits(j).V_L = data(:,k);
            bem.subunits(j).V_R = data(:,k+1);
        end
        
    else
        if ~bem.silent
            fprintf('No data found. Exiting.\n');
        end
        
    end
   
end