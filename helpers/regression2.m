function [r,m,b] =  regression2(y,x,noOffset)
    % Uage: [r,m,b] = regression2(y,x);
    % r is the r statistic
    % m is the slope
    % b is the offset
    % Fits a line : y = b + m*x, using least squares
    
    if nargin < 3
        noOffset = 0;
    end
    
    % Just try to sort the inputs out that so you get column vectors
    x = ascolumn(x);
    y = ascolumn(y);
        
    if any(isnan(y(:)));
        r = NaN;
        m = NaN;
        b = NaN;
        return
    end
    
    if ~noOffset        
        X = [ones(size(x,1),1),x];
    else
        X = [zeros(size(x,1),1),x];
    end
                
    betaHat = pinv(X' * X) * X' * y;
        
    b = betaHat(1);
    
    m = betaHat(2:end);
    
    R = corrcoef(X*betaHat, y);
    
    r = R(1,2);
    
    
end