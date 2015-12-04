function G = gabor(x,y,rf)  
    % Compute two-dimensional vertically oriented Gabor
    % Usage: G = gabor(x,y,rf);
    % x and y: linear arrays; for 2D RFs use meshgrid.
    % rf: RF structure; use rf=bem.default_rf_params() for template
    
    
    G = spatial_gaussian(x,y,rf) .* cos(2*pi*rf.f*(x-rf.x0) + rf.phi);
end