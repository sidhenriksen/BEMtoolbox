function twopassData = get_twopass_data(NimStruct)
    % Wrapper for running the two-pass analysis on a NimStruct.
    % This function returns a data structure with 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Constant parameters %%%
    
    RUNDUR = 30;
        
    LAG = 3;
    
    WINDOWSIZE = 4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    [Vint,Vtotal,tc] = run_2pass(NimStruct,RUNDUR,WINDOWSIZE,LAG);


    twopassData.internalVariance = Vint;
    
    twopassData.totalVariance = Vtotal;
    
    twopassData.twopassSpikeCount = tc;
        
end