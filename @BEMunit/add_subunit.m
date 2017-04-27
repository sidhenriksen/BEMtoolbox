function varargout = add_subunit(bem,varargin)
    % Add a new binocular subunit to the model
    % Usage: bem = bem.add_subunit(<params>)
    % Example:
    % bem = bem.add_subunit('L_phi',pi/2, 'R_phi',0);
    % Will add a subunit with left and right eye's phases of pi/2 and 0,
    % respectively. 

    rf = BEMunit.default_rf_params();

    % Loop over and set the field to the appropriate value if the
    % field is valid; throw an assertion error otherwise.
    for j = 1:2:(length(varargin)-1);
        if strcmp('L_',varargin{j}(1:2))
            assert(isfield(rf.left,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            rf.left.(varargin{j}(3:end)) = varargin{j+1};
        elseif strcmp('R_',varargin{j}(1:2));
            assert(isfield(rf.right,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            rf.right.(varargin{j}(3:end)) = varargin{j+1};
        else
            if ~isfield(rf,varargin{j});
                assert(isfield(rf.left,varargin{j}),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            end
            rf.(varargin{j}) = varargin{j+1};
        end
    end

    
    
    bem.n_subunits = bem.n_subunits+1;
    
    if isempty(fieldnames(bem.subunits));
        sub_k = 1;
    else
        sub_k = length(bem.subunits)+1;
        bem.subunits(sub_k) = bem.subunits(1);
    end

    % This will handle issues with the receptive field centre not being
    % set. If neither bem.x0 or subunit.x0 is set, then it will
    % use existing centres if available and throw a warning. If no
    % subunits exist, then it will generate a random one and throw
    % an error.
    if isempty(bem.x0) && ~isfield(rf.left,'x0');
        warned = 0;
        if ~isempty(fieldnames(bem.subunits));
            if isfield(bem.subunits(1).rf_params.left,'x0')
                if ~bem.silent
                    warning('No RF centre specified; using centre from first subunit');
                end
                warned = 1;
            else
                error('No RF centres found. This is not supposed to happen.')
            end

            if isfield(bem.subunits(1).rf_params.left,'y0')
                if ~warned && ~bem.silent;
                    warning('No RF centre specified; using centre from first subunit');
                end
            else
                error('No RF centres found. This is not supposed to happen.')
            end

        else
            bem.subunits(1) = struct();
            bem = bem.set_centers();
            if ~bem.silent
                warning('No RF centres specified; using random centres');
            end
        end

        x0_L = bem.subunits(1).rf_params.left.x0;
        y0_L = bem.subunits(1).rf_params.left.y0;
        x0_R = bem.subunits(1).rf_params.right.x0;
        y0_R = bem.subunits(1).rf_params.right.y0;

    else
        bem = bem.set_centers();

        x0_L = bem.subunits(sub_k).rf_params.left.x0;
        y0_L = bem.subunits(sub_k).rf_params.left.y0;

        x0_R = bem.subunits(sub_k).rf_params.right.x0;
        y0_R = bem.subunits(sub_k).rf_params.right.y0;
    end



    % Create arrays/grid for computing spatial RF
    x_max = bem.Nx * bem.deg_per_pixel;
    y_max = bem.Ny * bem.deg_per_pixel;            
    bem.x = linspace(0,x_max,bem.Nx);
    bem.y = linspace(0,y_max,bem.Ny);            
    [X,Y] = meshgrid(bem.x,bem.y);


    phi_L = -bem.dphi/2;
    phi_R = bem.dphi/2;

    rf.left.phi = rf.left.phi + phi_L;
    rf.left.x0 = x0_L;
    rf.left.y0 = y0_L;

    rf.right.phi = rf.right.phi + phi_R;
    rf.right.x0 = x0_R;
    rf.right.y0 = y0_R;


    % Create left eye receptive field
    switch rf.left.rf_type
        case 'gabor'
            bem.subunits(sub_k).L = gabor(X,Y,rf.left);
        case 'gaussian'
            bem.subunits(sub_k).L = spatial_gaussian(X,Y,rf.left);
    end

    % Create right eye receptive field
    switch rf.right.rf_type
        case 'gabor'
            bem.subunits(sub_k).R = gabor(X,Y,rf.right);
        case 'gaussian'
            bem.subunits(sub_k).R = spatial_gaussian(X,Y,rf.right);
    end

    bem.t = linspace(0,1.0,1001);
    tc = median(bem.t);

    % Left temporal kernel
    switch rf.left.temporal_kernel;
        case 'gamma_cosine'
            bem.subunits(sub_k).L_tk = gamma_cosine(bem.t,tc,rf.left);
        case 'gaussian'
            bem.subunits(sub_k).L_tk = gaussian(bem.t,tc,rf.left);
        case 'none'
            bem.t=[];
    end

    % Right temporal kernel
    switch rf.right.temporal_kernel;
        case 'gamma_cosine'
            bem.subunits(sub_k).R_tk = gamma_cosine(bem.t,tc,rf.right);
        case 'gaussian'
            bem.subunits(sub_k).R_tk = gaussian(bem.t,tc,rf.right);
        case 'none'
            bem.t=[];
    end

    bem.subunits(sub_k).NL = @(x)(x.^2);
    
    bem.subunits(sub_k).weight = 1;

    bem.subunits(sub_k).rf_params = rf;
    
    bem.update();
    
    if nargout
        varargout = {bem};
    end
end