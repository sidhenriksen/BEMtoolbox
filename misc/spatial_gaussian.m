function G = spatial_gaussian(x,y,rf)
    % G = spatial_gaussian(x,y,rf);
    % Compute two-dimensional Gabor
    % x and y: linear arrays, use meshgrid
    % rf : default rf parameters
        
    G = exp(-((x-rf.x0).^2 ./ (2*rf.sx^2) + ((y-rf.y0).^2)./(2*rf.sy^2)));
end
