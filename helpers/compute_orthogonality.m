function orthogonality = compute_orthogonality(X)
    % Returns a metric of orthogonality for a matrix X. If all columns
    % of X are orthogonal, this function will return 1. If only a single
    % orthogonal component is present, the function will return 0 (i.e.
    % if the rank is 1).
    % 
    % Parameters
    % -----------
    % X : N x K matrix; orthogonality is computed across the K columns
    %
    % Returns
    % --------
    % orthogonality : real number between 0 and 1. 

    if size(X,2) == 1
        orthogonality = 1;
        return
    end

    [~,~,~,~,v] = pca(X);
    
    N = size(X,2);
        
    AUC = sum(cumsum(v))/(N*100);
    
    orthogonality = 1 - (2*AUC -1);
    
    orthogonality = orthogonality * 1/(1-1/length(v));

end