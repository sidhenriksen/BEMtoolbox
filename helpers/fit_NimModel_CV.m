function NimModel = fit_NimModel_CV(NimCell,nExc,nInh,optStruct)
    % runs fit_NimModel with 4-fold cross-validation
    
    if nargin< 4
        optStruct = struct();
    end
    
    K = 5;
    
    [trainIndices,cvIndices] = kfold_split(NimCell,K);
    
    rEst = zeros(1,length(NimCell.dxs));
    for j = 1:K
        
        optStruct.indexTrain = trainIndices{j};
        optStruct.indexCv = cvIndices{j};
        
        NimModel = fit_NimModel(NimCell,nExc,nInh,optStruct);
        
        rEst(cvIndices{j}) = NimModel.rEst(cvIndices{j});
        
    end
    
    
    NimModel.rEst = rEst;
    

end

function [trainIndices,cvIndices] = kfold_split(NimCell,K)


    trials = cumsum([0,diff(NimCell.times)<0])+1;
    
    uniqueTrials = shuffle(unique(trials));

    % so now we split uniqueTrials into five equal bits
    
    trainIndices = cell(1,K);
    cvIndices = cell(1,K);
    
    nInFold = floor(length(uniqueTrials)/K);
    
    for k = 1:K;
        currentStart = 1 + (k-1)*nInFold;
        
        
        if k == K
            currentStop = length(uniqueTrials);
            
        else
            currentStop = k*nInFold;
                   
        end
        
        currentIndices = arrayfun(@(k)find(trials==k),uniqueTrials(currentStart:currentStop),...
            'uniformoutput',false);
        
        currentIndices = cell2mat(currentIndices);
        
        cvIndices{k} = ascolumn(currentIndices(:));
        trainIndices{k} = ascolumn(setdiff(1:length(trials),currentIndices(:)));
        
    end
    
    
end