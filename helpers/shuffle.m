function y = shuffle(x,dim)
    % Shuffles input array x along dimension dim
    % y = shuffle(x)
    if nargin < 2
        dim = 1;
    end
        
    x = squeeze(x);    
    
    if length(size(x)) < dim; % in case we call shuffle(x,5) on an array; do nothing.
        
        y = x;
        
    elseif any(size(x)==1) % then this is an array
        
        y = x(randperm(length(x)));
        
    else
        k = size(x,dim);
        idx = randperm(k);
        if dim == 1;
            y = x(idx,:,:,:);
        elseif dim ==2
            y = x(:,idx,:,:);
        elseif dim == 3;
            y = x(:,:,idx,:);
        elseif dim ==4;
            y = x(:,:,:,idx);
        else
            error('Unsupported dimension.');
        end
            
    end
        
end