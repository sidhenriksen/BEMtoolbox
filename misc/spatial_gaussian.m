function G = spatial_gaussian(x,y,rf)
    
    % Compute two-dimensional Gabor
    % x and y: linear arrays, use meshgrid
    % x0 and y0: horizontal and vertical centres, respectively.
    % sx and sy: horizontal and vertical SDs of Gaussian
    G = exp(-((x-rf.x0).^2 ./ (2*rf.sx^2) + ((y-rf.y0).^2)./(2*rf.sy^2)));
end
