function twopassData = get_twopass_data(NimStruct)
    % Wrapper for running the two-pass analysis on a NimStruct.
    % This function returns a data structure with 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Constant parameters %%%
    
    RUNDUR = 30;
        
    LAG = 4;
    
    WINDOWSIZE = 3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    [Vtotal,Vint,tc,sqrtMeans,sqrtVars] = run_2pass(NimStruct,RUNDUR,WINDOWSIZE,LAG);


    twopassData.internalVariance = Vint;
    
    twopassData.totalVariance = Vtotal;
    
    twopassData.twopassSpikeCount = tc;
    
    twopassData.twopassSqrtSpikeCount = sqrtMeans;
    
    twopassData.totalSqrtVariance = sqrtVars;
    
end