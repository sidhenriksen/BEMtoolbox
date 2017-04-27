function G = gabor(x,y,rf)  
    % Compute two-dimensional vertically oriented Gabor
    % Usage: G = gabor(x,y,rf);
    % x and y: linear arrays; for 2D RFs use meshgrid.
    % rf: RF structure; use rf=bem.default_rf_params() for template
    
    if nargin < 3
        rfboth = BEMunit.default_rf_params();
        rf = rfboth.left;
    end
    
    xp = x*cos(rf.theta) + y*sin(rf.theta);
    yp = -x*sin(rf.theta) + y*cos(rf.theta);
    
    %xp = (x-rf.x0)*cos(rf.theta) + (y-rf.y0)*sin(rf.theta);
    %yp = -(x-rf.x0)*sin(rf.theta) + (y-rf.x0)*cos(rf.theta);
        
    % override standard gaussian settings since xp and yp have already
    % been translated
    x0 = rf.x0*cos(rf.theta) + rf.y0 * sin(rf.theta);
    y0 = -rf.x0*sin(rf.theta) + rf.y0 * cos(rf.theta);
    rf.x0 = x0; rf.y0 = y0;
    
    G = spatial_gaussian(xp,yp,rf) .* cos(2*pi*rf.f*(xp-rf.x0) + rf.phi);
end