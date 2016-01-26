function save_bootstrap(bem,generator)
    % Saves the computed monocular responses (pre-temporal
    % filtering) to a file with a unique identifier for this model
    % + stimulus. Will extend the file if it already exists.,
    % Usage: bem.save_bootstrap(generator);
    % generator: a stimulus generator object

    bem_id = bem.get_identifier();
    stim_id = generator.get_identifier();
    
    csvfile = [bem.bootstrap_dir,num2str(bem_id),'_',num2str(stim_id),'.csv'];
        
    new_data = zeros(length(bem.subunits(1).V_L),length(bem.subunits)*2);
        
    for j = 1:length(bem.subunits);
        
        k=(j-1)*2+1;
        new_data(:,k:k+1) = [bem.subunits(j).V_L,bem.subunits(j).V_R];
    end
    
    if exist(csvfile,'file');
        
        existing_data=csvread(csvfile);
        
        new_data = [existing_data;new_data];
    end
    
    
    csvwrite(csvfile,new_data);
    
end