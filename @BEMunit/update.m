function bem = update(bem)
    % Updates the bem_unit after changes have been made to RF
    % params

    % Loop over and set the field to the appropriate value if the
    % field is valid; throw an assertion error otherwise.

    % Create arrays/grid for computing spatial RF
    x_max = bem.Nx * bem.deg_per_pixel;
    y_max = bem.Ny * bem.deg_per_pixel;            
    bem.x = linspace(-x_max/2,x_max/2,bem.Nx);            
    bem.y = linspace(-y_max/2,y_max/2,bem.Ny);            
    [X,Y] = meshgrid(bem.x,bem.y);


    for j = 1:length(bem.subunits);

        % This allows you to set the bem.dx,bem.dy and bem.dphi
        % fields to be empty ([]) and use custom ones for each
        % subunit.
        if isempty(bem.dx)
            dx_L = bem.subunits(j).dx/2;
            dx_R = -bem.subunits(j).dx/2;
        else
            dx_L = bem.dx/2;
            dx_R = -bem.dx/2;
        end

        if isempty(bem.dy);
            dy_L = bem.subunits(j).dy/2;
            dy_R = -bem.subunits(j).dy/2;
        else
            dy_L = bem.dy/2;
            dy_R = -bem.dy/2;
        end

        if isempty(bem.dphi)
            phi_L = bem.subunits(j).dphi/2;
            phi_R = -bem.subunits(j).dphi/2;
        else
            phi_L = bem.dphi/2;
            phi_R = -bem.dphi/2;
        end




        rf = bem.subunits(j).rf_params;

        if ~isempty(bem.x0)
            rf.left.x0 = bem.x0 + dx_L;
            rf.right.x0 = bem.x0 + dx_R;
        end

        if ~isempty(bem.x0)
            rf.left.y0 = bem.y0 + dy_L;
            rf.right.y0 = bem.y0 + dy_R;
        end



        rf.left.phi = rf.left.phi + phi_L;
        rf.right.phi = rf.right.phi + phi_R;

        % Create left eye receptive field
        switch rf.left.rf_type
            case 'gabor'
                bem.subunits(j).L = gabor(X,Y,rf.left);
            case 'gaussian'
                bem.subunits(j).L = spatial_gaussian(X,Y,rf.left);
        end

        % Create right eye receptive field
        switch rf.right.rf_type
            case 'gabor'
                bem.subunits(j).R = gabor(X,Y,rf.right);
            case 'gaussian'
                bem.subunits(j).R = spatial_gaussian(X,Y,rf.right);
        end

        bem.t = linspace(0,1.0,1001);
        t_c = median(bem.t);

        % Left temporal kernel
        switch rf.left.temporal_kernel;
            case 'gamma_cosine'
                bem.subunits(j).L_tk = gamma_cosine(bem.t,t_c,rf.left);
            case 'gaussian'
                bem.subunits(j).L_tk = gaussian(bem.t,t_c,rf.left);
            case 'none'
                bem.t=[];
        end

        % Right temporal kernel
        switch rf.right.temporal_kernel;
            case 'gamma_cosine'
                bem.subunits(j).R_tk = gamma_cosine(bem.t,t_c,rf.right);
            case 'gaussian'
                bem.subunits(j).R_tk = gaussian(bem.t,t_c,rf.right);
            case 'none'
                bem.t=[];
        end
        bem.subunits(j).rf_params = rf;

    end
    bem.n_subunits = length(bem.subunits);
end