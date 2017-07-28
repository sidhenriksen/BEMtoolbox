function [r,m,b] = type2_regression(x,y,lambda)    
     % type2_regression
     % Wrapper for Jenny's fit_bothsubj2error function (type 2 regression
     % function) in order to make
     % output consistent with Matlab's regression method.
     % Usage:
     % [r,m,b] = type2_regression(x,y,lambda)
     %
     % Parameters
     % -----------
     % x,y: input arrays
     %  
     % lambda : scalar
     %      the ratio of the variances for y and x
     %      lambda = sum(var[y]) / sum(var[x])
     %
     % Returns
     % ---------
     % r : correlation coefficient
     % m : slope 
     % b : offset
     
    if nargin <3
        lambda = Inf;
    end
    
    [b,m] = fit_bothsubj2error(x,y,lambda);
    
    r = corrcoef(b + m*x,y);
    r=r(2);

end