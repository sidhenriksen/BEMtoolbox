function [r,m,b] =  regression2(y,x)
    % Uage: [r,m,b] = regression2(y,x);
    % r is the r statistic
    % m is the slope
    % b is the offset
    % Fits a line : y = b + m*x, using least squares
    
    % Just try to sort the inputs out that so you get column vectors
    if ~iscolumn(x)
        x = x';
    end
    if ~iscolumn(y)        
        y = y';
    end        
    
    if any(isnan(y(:)));
        r = NaN;
        m = NaN;
        b = NaN;
        return
    end
    
    X = [ones(size(x,1),1),x];
                
    betaHat = pinv([X' * X]) * X' * y;
    
    b = betaHat(1);
    
    m = betaHat(2:end);
    
    R = corrcoef(X*betaHat, y);
    
    r = R(1,2);
    
    
end