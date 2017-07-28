function bestParams = genetic_kernel_opt(NimModel,initParams,optParams)
    % GENETIC_KERNEL_OPT    Genetic kernel optimisation for separating high
    % vs low contrast responses.
    % 
    % Parameters
    % -----------
    % NimModel : NimModel structure
    % initParams : initial parameters
    %       [aL,mL,sL,fL,pL,aR,mR,sR,fR,pR], where a,m,s,f, and p
    %       refer to amplitude, mean, SD, freq, and phase, respectively.
    %       L and R refer to left vs right. Omit f and p for Gaussian.
    % optParams : a struct of optional parameters
    % (silent,fertilityRate,mutationRate,nGenerations)
    %
    % Returns:
    % bestParams : best parameter vector found
    
    silent = 0;
    
    fertilityRate = 20; % how many offspring per generation, integer

    mutationRate = 0.125; % how different the offspring are from the parent
    
    nGenerations = 30; % number of generations
    
    elitistSelection = true;
    
    costType = 2; 
    
    plotProgress = 0; 
    
    
    
    %pTrials = 0.8; % proportion training data

    modelType = 'gabor';
    
        
    if nargin > 2
        unpackStruct(optParams);
    end

    if strcmp(modelType,'gaussian')
        initParams([4:5,9:10]) = 0;                
        
    end
        
    bestParams = initParams;
    
    parentParams = initParams;
    
    
    bestCost = 0;
    
    if ~silent
        fprintf('\n Running genetic search\n')
    end
    
    if plotProgress;
        myFig = figure();
    end
    
    if ~isfield(NimModel,'seed');
        NimModel.seed = randi(1e7);
    end
    
    for j = 1:nGenerations;
        
        if ~silent
            fprintf('Generation %i of %i\n',j,nGenerations);
        end
        
        [candidateBestParams,candidateBestCost] = find_best_genetic_params(NimModel,parentParams,mutationRate,fertilityRate,modelType,costType);
                
        if candidateBestCost > bestCost
                        
            bestCost = candidateBestCost;
            
            bestParams = candidateBestParams;
                        
            parentParams = bestParams;
                   
        end
        
        if plotProgress
            figure(myFig); cla;
            [thisCost,bestKernel] = evaluate_params(NimModel,bestParams,0,costType);
            plot(bestKernel)          
            
            title(sprintf('Generation %i of %i; cost=%.4f',j,nGenerations,thisCost));
            drawnow;
        end
        
        

        % excluding this means we have elitist selection
        if ~elitistSelection;
            parentParams = candidateBestParams;
        end
                
    end
    
end


function [bestParams,bestCost] = find_best_genetic_params(NimModel,parentParams,mutationRate,fertilityRate,modelType,costType)

    is2d = strcmpi(modelType,'gaussian2d');
    
    rng('shuffle');

    mutationParams = parentParams; 
    
    if ~all(parentParams([4:5,9:10]) == 0);
        mutationParams([5,10]) = pi; % set phase to pi if 
    end
    
    if is2d
        monocularMutationParams = [0.15,5,5,2,2,0.25];
        mutationParams = [monocularMutationParams,monocularMutationParams];
    end

    mutationParams([1,6]) = 2;
    
    currentMutationRates = repmat(mutationParams*mutationRate,[fertilityRate,1]);
    

        
    childParams = repmat(parentParams,[fertilityRate,1]) + ...
            randn(fertilityRate,length(parentParams)).*currentMutationRates;
        
    if is2d % this just sets values below -1 or above 1 to -1 and 1, respectively
        leftRect = abs(childParams(:,6)) > 1;
        leftSign = sign(childParams(leftRect,6));
        childParams(leftRect,6) = leftSign;
        
        rightRect = abs(childParams(:,12)) > 1;
        rightSign = sign(childParams(rightRect,12));
        childParams(rightRect,12) = rightSign;
    end
        
        
    costs = arrayfun(@(k)evaluate_params(NimModel,childParams(k,:),is2d,costType),1:fertilityRate);
    
    [bestCost,bestParamIdx] = max(costs);
    
    bestParams = childParams(bestParamIdx,:);
    %fprintf('Best cost: %.3f\n',bestCost);

end