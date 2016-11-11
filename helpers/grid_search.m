function [nExc,nInh] = grid_search(NimCell)
   % Do exhaustive grid search on many combinations of inhibitory and
   % excitatory subunits.
   % Usage:
   % [nExc,nInh] = grid_search(NimCell);
   % NimCell is a NimCell object (duh)
   % nExc and nInh are the number of excitatory and inhibitory subunits,
   % respectively.

    nRepeats = 3;
    
    nExcs = 1:7;
    
    nInh = 0:5;        

    paramGrid = CombVec(nExcs,nInh)';
    
    optStruct.pTrain = 0.75; % proportion of data used for training set    
    optStruct.silent = 1;
        
    LLs = zeros(size(paramGrid,1),nRepeats); % log likelihoods    
        
    for k = 1:size(paramGrid,1);
        
        params = paramGrid(k,:);
        
        nExc = params(1);
        
        nInh = params(2);
        
        fprintf('Running params: nExh=%i, nInh=%i\n',nExc,nInh);
        for j = 1:nRepeats;
            
            % returns cross-validated log likelihood            
            [~,LLCv] = fit_NimModel(NimCell,nExc,nInh,optStruct);
            
            LLs(k,j) = LLCv;
            
        end
                
    end
    
    mLLs = median(LLs,2); % median or mean
    
    [~,idx] = max(mLLs);
    
    nExc = paramGrid(idx,1);
    
    nInh = paramGrid(idx,2);

end