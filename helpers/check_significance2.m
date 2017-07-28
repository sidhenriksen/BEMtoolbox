function P = check_significance2(NimCell,N)

    if nargin < 3
        N = 1e5;
    end
    
    runDur = 30;
    windowSize = 1;
    lag0 = 0;
    lag1 = get_maxvar_time(NimCell);
    

    [~,tc0] = run_forward_correlation(NimCell,runDur,windowSize,lag0);
    [~,tc1] = run_forward_correlation(NimCell,runDur,windowSize,lag1);
    
    cTc0 = tc0(3,:);
    cTc1 = tc1(3,:);
    
    P = permutation_test(cTc1,cTc0,N);
end


function h = permutation_test(x,y,N)

    if nargin < 3
        p = 0.01;
    end
    
    if nargin < 4
        N = 1e4;
    end
    
    a = x-mean(x);
    b = y-mean(y);
    permutationDifferences = arrayfun(@(k)permutation_difference(a,b),1:N);
    
    realDifference = var(x) - var(y);
    
    h = mean( realDifference <= permutationDifferences);
    
end

function d = permutation_difference(x,y)

    z = shuffle([x,y]);
    d = var(z(1:length(x))) - var(z(length(x)+1:end));
    
end