function G = gabor(x,y,rf)  
    % Compute two-dimensional vertically oriented Gabor
    % Usage: G = gabor(x,y,rf);
    % x and y: linear arrays; for 2D RFs use meshgrid.
    % rf: RF structure; use rf=bem.default_rf_params() for template
    
    xp = x*cos(rf.theta) + y*sin(rf.theta);
    yp = -x*sin(rf.theta) + y*cos(rf.theta);
    
    G = spatial_gaussian(xp,yp,rf) .* cos(2*pi*rf.f*(xp-rf.x0) + rf.phi);
end