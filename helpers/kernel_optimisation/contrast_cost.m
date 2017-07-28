function cost = contrast_cost(NimModel,binocularKernel,costType)
    % CONTRAST_COST Evaluates the current binocular kernel,
    % and returns the cost in mean spike count difference
    % for high and low contrast frames.
    %
    % Parameters
    % -----------
    % NimModel : NimModel structure
    % binocularKernel : a 1-D array whose length is size(NimModel.stim,2)
    %
    % Returns
    % --------
    % cost : mean spike count difference for high vs low contrast frames.
 
    if isfield(NimModel,'seed')
        optStruct.seed = NimModel.seed;
    else
        
        rng('shuffle')
        
        optStruct.seed = randi(1e6);
        
    end
    
    if isfield(NimModel,'cv');
        
        optStruct.cv = NimModel.cv;
                
    else
        optStruct.cv = 0;
    end
    
    
    if nargin < 3
        costType = 2;
    end
           
    optStruct.pSplit = 1/5;
    
    optStruct.runDur = 30;
                
    
    [lo,hi] = split_trials_by_contrast(NimModel,binocularKernel,optStruct);
    
    if ~sum(lo) || ~sum(hi)
        
        cost = 0; return;
        
    end
    
    step = 3;
    
    if length(binocularKernel)>500
        lag = 1;
    else
        lag = 4;
    end
    
    if isfield(NimModel,'rObs')
        rObs = NimModel.rObs;
    elseif isfield(NimModel,'rEst')
        rObs = NimModel.rEst;    
    else
        warning('No rates on structure. This gets computed locally, which adds a fair amount of overhead. Consider pre-computing.')
        rObs = get_rate(NimModel);
    end
    
    
    if costType == 1
        dxs = unique(NimModel.dxs);
        dxs = dxs(abs(dxs) < 10);
        corrs = [-1,0,1];



        rLo = zeros(length(corrs),length(dxs));

        rHi = zeros(length(corrs),length(dxs));



        for c = 1:length(corrs)
            currentCorrIdx = (corrs(c) == NimModel.correlation);
            for k = 1:length(dxs)

                if c == 2;
                    currentDxIdx = currentCorrIdx;

                else
                    currentDxIdx = (dxs(k) == NimModel.dxs);
                end

                currentLoIdx = lo & currentCorrIdx & currentDxIdx;
                currentHiIdx = hi & currentCorrIdx & currentDxIdx;

                currentRLo = get_spike_counts(rObs,currentLoIdx,lag,step);

                currentRHi = get_spike_counts(rObs,currentHiIdx,lag,step);



                if c == 2;
                    rLo(c,:) = mean(currentRLo);
                    rHi(c,:) = mean(currentRHi);
                    break;
                end

                rLo(c,k) = mean(currentRLo);
                rHi(c,k) = mean(currentRHi);

            end
        end
        
        [~,m,~] = regression2(rHi(:)',rLo(:)');

        %cost = mean(rHi)-mean(rLo);
    
        cost = m;

    elseif costType == 2
        
       rLo = get_spike_counts(rObs,lo,lag,step);
       
       rHi = get_spike_counts(rObs,hi,lag,step);
       
       cost = mean(rHi)-mean(rLo);
       
    end
    
end