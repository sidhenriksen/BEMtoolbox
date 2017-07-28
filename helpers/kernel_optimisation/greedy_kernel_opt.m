function bestParams = greedy_kernel_opt(NimModel,initParams,threshold,optParams)

    silent = 0;
    nMaxHard = 100;
    nMaxSoft = 20;
    modelType = 'gabor';

    if nargin > 3
        unpackStruct(optParams);
    end

    switch modelType
        case 'gabor'
            paramsToSearch = {'a_L','m_L','s_L','f_L','p_L','a_R','m_R','s_R','f_R','p_R'};
        case 'gaussian'
            paramsToSearch = {'a_L','m_L','s_L','a_R','m_R','s_R'};
            initParams([4:5,9:10]) = 0;
            
        case 'gaussian2d'
            paramsToSearch = {...
                'a_xL','m_xL','s_xL','m_tL','s_tL','r_L'...
                'a_xR','m_xR','s_xR','m_tR','s_tR','r_R'};
                
    end
    
    %assert(length(initParams)==length(paramsToSearch),'Invalid parameter length.');
    
    k = 0;
    
    bestParams = initParams;
    
    bestCost = 0;
    
    thresholdFailureCounter = 0;
    
    if ~silent
        fprintf('Running greedy search.\n');
    end
    
    for j = 1:nMaxHard;
        k = k+1;
        
        if ~mod(k,length(paramsToSearch)+1)
            k = 1;
        end
        
                        
        paramToSearch = paramsToSearch{k};
        
        [candidateBestParams,candidateBestCost] = find_best_greedy_params(NimModel,bestParams,paramToSearch,modelType);
         
        if candidateBestCost < (bestCost+threshold) % if no significant improvement
            
            thresholdFailureCounter = thresholdFailureCounter +1;
            
            if ~silent
                fprintf('Failed to improve; %i of %i\n',thresholdFailureCounter,nMaxSoft);
            end
            
        else
            
            bestCost = candidateBestCost;
            
            bestParams = candidateBestParams;
            
        end
        
        if thresholdFailureCounter >= nMaxSoft
            
            break
            
        end                
    end
        
end

function [bestParams,bestCost] = find_best_greedy_params(NimModel,currentParams,binocularParamToSearch,modelType)

    if nargin < 4
        modelType = 'gabor';
    end

        
    searchList = {...
        linspace(0,2,21),... % amplitudes
        [1:10,10:3:20,20:30,30:3:51], ... % means
        linspace(0.5,10,21), ... % SDs
        0, ...%linspace(0.01,0.2,21), ... % freqs
        0, ...%linspace(0,2*pi,11) ... % phases    
        linspace(-1,1,21) ... % correlation (gaussian tilt)
    };

    is2d = strcmpi(modelType,'gaussian2d');

    paramList = {'a','m','s','f','p','r'};
    
    whichEye = binocularParamToSearch(end);
    
    paramToSearch = binocularParamToSearch(1);
    
    paramIdx = strcmp(paramList,paramToSearch);
    
    binocularParamIdx = find(paramIdx) + strcmp(whichEye,'R')*5; 
    
    currentSearchArray = searchList{paramIdx};
    
    paramCombinations = repmat(currentParams,[length(currentSearchArray),1]);
    
    paramCombinations(:,binocularParamIdx) = currentSearchArray;
    
    % get cost for each of the parameter combinations
    costs = arrayfun(@(k)evaluate_params(NimModel,paramCombinations(k,:),is2d),1:length(currentSearchArray));
    
    [bestCost,bestParamIdx] = max(costs);
    
    bestParams = paramCombinations(bestParamIdx,:);
    
    
end
