function NimModel = generate_new_stimuli(NimModel)
    % NimModel = generate_new_stimuli(NimModel)
    % Takes the NimModel and generates a new sequence of 1D RDS by
    % resampling from existing stimuli
    

    N = size(NimModel.stim,1);
    K = round(N/3);
    
    replacement=1;
    
    idx = randsample(N,K,replacement);
    idx = repmat(idx,[1,3])';
    idx = idx(:);
            
    NimModel.stim = NimModel.stim(idx,:);
    NimModel.correlation = NimModel.correlation(idx);
    NimModel.duration = NimModel.duration*0 + 30;
    NimModel.dxs = NimModel.dxs(idx);
           

end