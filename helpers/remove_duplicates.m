function Y = remove_duplicates(X,k)
    % Function to remove duplicate numbers; if k=3, and sequence is
    % S = [i,i+1,i+2], remove_duplicates will return an array with
    % i+1 and i+2 removed.
    % Usage: Y = remove_duplicates(X,k);

    skipCounter = k-1;
    
    k = round(k);   
    
    keepIndex = zeros(size(X));
    for j = 1:length(X);
        
        if skipCounter < (k-1)
            skipCounter = skipCounter+1;
            continue;
        end
        
        
        stop = min(j+k-1,length(X));
        if sum(diff(X(j:stop))) == k-1;
            keepIndex(j)=1;
            skipCounter = 0;
        end
        
    end
    
    Y = X(logical(keepIndex));
    
end
