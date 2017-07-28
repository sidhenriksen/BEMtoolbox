function contrastData = get_contrast_data(NimStruct,binocularKernel)

    if nargin < 2
        % or just set to empty
        binocularKernel = kernel_optimisation(NimStruct,'genetic');
                
    end


    assert(isfield(NimStruct,'NimFit'),'get_contrast_data needs NimFit field to run')


    %%% Modifiable parameters %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    optParams.pSplit = 1/3;    
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Constant parameters %%%
    
    RUNDUR = 30;
    
    WINDOWSIZE = 3;
    
    LAG = 4;
        
    triggerOnce = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    [lo,hi] = split_trials_by_contrast(NimStruct,binocularKernel,optParams);

    % this gets the tuning curves and variances using forward correlation
    [dx,tcLo,varLo] = run_forward_correlation(NimStruct,RUNDUR,WINDOWSIZE,LAG,triggerOnce,lo);

    [~,tcHi,varHi] = run_forward_correlation(NimStruct,RUNDUR,WINDOWSIZE,LAG,triggerOnce,hi);

    
    % Stick the data into a data structure, and then return the structure
    
    contrastData.dx = dx;
    
    contrastData.spikeCountLowContrast = tcLo;
    
    contrastData.spikeCountHighContrast = tcHi;
    
    contrastData.varianceLowContrast = varLo;
    
    contrastData.varianceHighContrast = varHi;
    
    contrastData.contrastKernel = binocularKernel;



end