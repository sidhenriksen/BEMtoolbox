% Class for creating model binocular cells
% Author: Sid Henriksen
% Laboratory of Sensorimotor Research - National Institutes of Health (US)
% Institute of Neuroscience - Newcastle University (UK)
% Email: sid.henriksen@gmail.com
% 
% Please email me with any bugs or suggestions.
%
% All code released under GNU GPL v2
classdef BEMunit
    properties
        
        % spatial dimensions; valid is 1d and 2d
        dim='2d';
        
        x0;
        y0;
        
        % binocular properties
        dx=0; %horizontal disparity
        dy=0; %vertical disparity
        dphi=0; %phase disparity
        
        dt=0.001;
        
        %number of pixels
        Nx=200;
        Ny=200;
        
        deg_per_pixel=0.025;
        
        subunits=struct();
        
        temporal_kernel='none';
        
    end
    
    properties (Hidden, SetAccess=private)
        x;
        y;
        t;
        n_subunits=2;
        
        init=true;
        
        bootstrap_dir;
        toolbox_dir;
        
    end
    
    properties (Hidden)
        memory_threshold=5e9; % in bytes; currently 5GB.
    end
    
    methods
        function bem = BEMunit(varargin)
            % Loop over and set the field to the appropriate value if the
            % field is valid; throw an assertion error otherwise.
            for j = 1:2:(length(varargin)-1);
                assert(isprop(bem,varargin{j}),'Error: Invalid field given to BEM constructor');
                bem.(varargin{j}) = varargin{j+1};
            end
        
            
            % it 'init' is set to false, don't initialise to a standard BEM
            % unit. Default is to initialise a standard BEM unit.
            if bem.init
                
                bem = add_subunit(bem);
                bem = add_subunit(bem,'L_phi',pi/2,'R_phi',pi/2);
                bem.toolbox_dir=mfilename('fullpath');
                bem.bootstrap_dir=[bem.toolbox_dir,'/.data/'];
            end
            
            
        end
        
        function C = simulate_spatial(bem,generator,n_frames)
            
            if nargin < 3;
                n_frames = 1;
            end
            
            
            S = zeros(bem.n_subunits,n_frames); % subunit responses
            
            total_bytes = bem.Nx * bem.Ny *n_frames * 8;
            if total_bytes > bem.memory_threshold;
                warning(['Dimensionality of stimulus too high; using serialised implementation to conserve memory (this is slower). ', ...
                    'You can manually increase the limit by changing bem.memory_threshold property.'])
                error('Serialisation not yet implemented. Reduce dimensionality or increase memory threshold');
                
            else
                        
                I_Ls = zeros(bem.Nx*bem.Ny,n_frames);
                I_Rs = zeros(bem.Nx*bem.Ny,n_frames);

                for j = 1:n_frames;
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(:);
                    I_Rs(:,j) = I_R(:);
                end    
                
                for k = 1:bem.n_subunits;
                    % compute left and right eye responses
                    L = bem.subunits(k).L(:)' * I_Ls;
                    R = bem.subunits(k).R(:)' * I_Rs;
                    
                    % add L and R, and pass through nonlinearity
                    S(k,:) = bem.subunits(k).NL(L+R);                     
                end
            end
            
            % Add together subunit responses to get final response
            % (complex cell response in the case of an energy model unit)
            C = sum(S);
            
        end
        
        function V = simulate_spatiotemporal(bem,generator,n_frames,duration,varargin)
            % Computes the spatiotemporal response of the cell in response
            % to a trial consisting of the frames given by the stimulus
            % generator. 
            save_bootstrap = false;
            load_bootstrap = false;
            for j = 1:length(varargin);
                switch varargin{j}
                    case 'save_bootstrap'
                        save_bootstrap = varargin{j+1};
                    case 'load_bootstrap'
                        load_bootstrap = varargin{j+1};
                end
            end
            
            if nargin < 3
                n_frames = 10;
            end
            
            if nargin < 4; 
                duration = 0.5;
            end
            
            
            n_frame_repeats = duration/(n_frames*bem.dt);           
            
            
            total_bytes = bem.Nx * bem.Ny *n_frames * 8;
            if total_bytes > bem.memory_threshold;
                warning(['Dimensionality of stimulus too high; using serialised implementation to conserve memory (this is slower). ', ...
                    'You can manually increase the limit by changing bem.memory_threshold property.'])
                error('Serialisation not yet implemented. Reduce dimensionality or increase memory threshold');
                
            else
                        
                I_Ls = zeros(bem.Nx*bem.Ny,n_frames);
                I_Rs = zeros(bem.Nx*bem.Ny,n_frames);

                for j = 1:n_frames;
                    [I_L,I_R] = generator.generate();
                    I_Ls(:,j) = I_L(:);
                    I_Rs(:,j) = I_R(:);
                end    
                
                Ls = zeros(bem.n_subunits,n_frames*n_frame_repeats);
                Rs = zeros(bem.n_subunits,n_frames*n_frame_repeats);
                
                
                for k = 1:bem.n_subunits;
                    % compute left and right eye responses
                    L = bem.subunits(k).L(:)' * I_Ls;
                    R = bem.subunits(k).R(:)' * I_Rs;
                    
                    L_rep = repmat(L,[1,n_frame_repeats]);
                    R_rep = repmat(R,[1,n_frame_repeats]);
                    
                    Ls(k,:) = L_rep;
                    Rs(k,:) = R_rep;
                end
            end
            
            if isempty(bem.bootstrap_dir);
                bem.bootstrap_dir = 
                
            end
            
            
            bem.t = linspace(0,(length(Vhat)-1)*bem.dt,length(Vhat));
            t_c = median(bem.t);
            switch bem.temporal_kernel
                case 'gaussian'
                    tk=temporal_gaussian(bem.t,t_c,bem.tau);
                    
                case 'gamma-cosine'
                    tk=gamma_cosine(bem.t,t_c,bem.alpha,bem.omega,bem.tau,bem.t_phi);
            end
            tk_fft = fft(tk);
            Vhat_fft = fft(Vhat);
            
            V = ifft(tk_fft .* Vhat_fft);
            
        end
        
        function bem = update(bem)
            % Updates the bem_unit after changes have been made to RF
            % params
            
            % Loop over and set the field to the appropriate value if the
            % field is valid; throw an assertion error otherwise.
            
            % Create arrays/grid for computing spatial RF
            x_max = bem.Nx * bem.deg_per_pixel;
            y_max = bem.Ny * bem.deg_per_pixel;            
            bem.x = linspace(-x_max/2,x_max/2,(bem.Nx));            
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
        
        function display_model(bem,plot_kernels)
            if nargin ==1
                plot_kernels = 0;
            end
            switch bem.dim
                case '2d'
                
                figure(); colormap('gray')
                
                % ticks for the axes
                xtk = linspace(min(bem.x),max(bem.x),5);
                ytk = linspace(min(bem.y),max(bem.y),5);
                clims = zeros(2*bem.n_subunits,2);
                
                for j = 1:bem.n_subunits;
                    subplot(bem.n_subunits,2,(j-1)*2+1);
                    imagesc(bem.x,bem.y,bem.subunits(j).L);
                    title(sprintf('Subunit %i: Left eye',j),'fontsize',14);
                    if j == 2;
                        xlabel('Horizontal position (deg)');
                    end
                    ylabel('Vertical position (deg)');
                    set(gca,'xtick',xtk,'ytick',ytk);
                    clims((j-1)*2+1,:)=get(gca,'clim');

                    subplot(bem.n_subunits,2,j*2)
                    imagesc(bem.x,bem.y,bem.subunits(j).R);
                    title(sprintf('Subunit %i: Right eye',j),'fontsize',14);
                    
                    if j == 2;
                        xlabel('Horizontal position (deg)');
                    end
                    
                    ylabel('Vertical position (deg)');
                    set(gca,'xtick',xtk,'ytick',ytk);
                end
                
                
                for j = 1:(bem.n_subunits*2);
                    subplot(bem.n_subunits,2,j);
                    set(gca,'clim',[min(clims(:)),max(clims(:))]);
                end
                
                if plot_kernels
                    figure();
                    subplot(2,1,1); hold on;
                    for j = 1:bem.n_subunits;
                        L = bem.subunits(j).L;
                        R = bem.subunits(j).R;
                        [~,max_idx] = max(var(L,[],1));

                        % cross section of left and right RF
                        xL = L(max_idx,:);
                        xR = R(max_idx,:);

                        % cross-correlation function
                        ccf=conv(xL,xR(end:-1:1));
                        cc_idx=round(length(xL)/2+1):round(length(xL)*1.5);
                        ccf = ccf(cc_idx);

                        xdx=cc_idx*bem.deg_per_pixel; xdx=xdx-median(xdx);

                        plot(xdx,ccf,'linewidth',3);
                        xlabel('Disparity (deg)','fontsize',14)
                        ylabel('Responses','fontsize',14);

                    end

                    subplot(2,1,2);
                    for j = 1:bem.n_subunits;
                        % Do something

                    end
                end
                case '1d'
                    % Plot 1d RF shere...
            end
        end
        
        function bem = add_subunit(bem,varargin)
            % Add a new binocular subunit to the model
            % Usage: 
            
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
                        warning('No RF centre specified; using centre from first subunit');
                        warned = 1;
                    else
                        error('No RF centres found. This is not supposed to happen.')
                    end
                    
                    if isfield(bem.subunits(1).rf_params.left,'y0')
                        if ~warned;
                            warning('No RF centre specified; using centre from first subunit');
                        end
                    else
                        error('No RF centres found. This is not supposed to happen.')
                    end
                    
                else
                    bem.subunits(1) = struct();
                    bem = bem.set_centers();
                    warning('No RF centres specified; using random centres');
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
            
            
            phi_L = bem.dphi/2;
            phi_R = -bem.dphi/2;
                
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
            
            bem.subunits(sub_k).rf_params = rf;
        end
        
        
        function bem = modify_subunit(bem,k_sub,varargin)
            % Usage: bem = bem.modify_subunit(k_sub,<attributes>)
            % Where k_sub is the index of the subunit to be modified
            %
            % ===Example usage===
            % Change left and right RF phase:
            % bem = bem.modify_subunit(2,'L_phi',pi/2,'R_phi',0);
            % 
            % Change output NL to half-squaring (here Pos would be a
            % pre-defined half-wave rectifying function):
            % bem = bem.modify_subunit(1,'NL',@(x)(Pos(x).^2));
            
            rf = bem.subunits(k_sub).rf_params;
            
            for j = 1:2:(length(varargin)-1);
                if strcmp('L_',varargin{j}(1:2))
                    assert(isfield(rf.left,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
                    rf.left.(varargin{j}(3:end)) = varargin{j+1};
                elseif strcmp('R_',varargin{j}(1:2));
                    assert(isfield(rf.right,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
                    rf.right.(varargin{j}(3:end)) = varargin{j+1};
                elseif strcmpi('NL',varargin{j});
                    rf.NL = varargin{j+1};                    
                else
                    if ~isfield(rf,varargin{j});
                        assert(isfield(rf.left,varargin{j}),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
                    end
                    rf.(varargin{j}) = varargin{j+1};
                end
            end
            
            bem.subunits(k_sub).rf_params = rf;
            bem = bem.update();
        end
        
        function bem = rescale(bem,scale)
            
            for j = 1:length(bem.subunits);
                
                bem.subunits(j).rf_params.left.sx=bem.subunits(j).rf_params.left.sx*scale;
                bem.subunits(j).rf_params.left.sy=bem.subunits(j).rf_params.left.sy*scale;
                bem.subunits(j).rf_params.left.f=bem.subunits(j).rf_params.left.f/scale;
                
                bem.subunits(j).rf_params.right.sx=bem.subunits(j).rf_params.right.sx*scale;
                bem.subunits(j).rf_params.right.sy=bem.subunits(j).rf_params.right.sy*scale;
                bem.subunits(j).rf_params.right.f=bem.subunits(j).rf_params.right.f/scale;
                
            end
            bem = bem.update();
        end

            
    end
    
   
    methods (Static)
       function rf_params = default_rf_params(~)
            % Usage: rf_params = default_rf_params();
            % Returns a struct with default rf parameters as fields

            %monocular (spatial) properties of receptive field
            rf_params_both.rf_type='gabor';
            %spatial extent of rf
            rf_params_both.sx=0.1;
            rf_params_both.sy=0.1;
            rf_params.sigma=[];
            %frequency of gabor
            rf_params_both.f=0.3./0.1;
            rf_params_both.phi=0;

            % temporal properties of receptive field
            rf_params_both.temporal_kernel='gamma-cosine';
            rf_params_both.tau=0.015;
            rf_params_both.omega=8.3;
            rf_params_both.t_phi=-pi/5;
            rf_params_both.alpha=2.5;
            

            % binocular properties
            rf_params.dphi=0; %phase disparity
            rf_params.dx=0; %position disparity (horizontal)
            rf_params.dy=0; %position disparity (vertical)

            rf_params.left=rf_params_both;
            rf_params.right=rf_params_both;
       end
    end
    
    methods (Access=private)
        function bem = set_centers(bem)
            x_max = bem.Nx * bem.deg_per_pixel;
            y_max = bem.Ny * bem.deg_per_pixel;
            
            if isempty(bem.x0);
                bem.x0 = (rand*0.3 + 0.3)*x_max;
            end
            if isempty(bem.y0);
                bem.y0 = (rand*0.3 + 0.3)*y_max;
            end

            x0_L = bem.x0 + bem.dx/2;
            x0_R = bem.x0 - bem.dx/2;

            y0_L = bem.y0 + bem.dy/2;
            y0_R = bem.y0 - bem.dy/2;
            
            for j = 1:length(bem.subunits);
                bem.subunits(j).rf_params.left.x0 = x0_L;
                bem.subunits(j).rf_params.right.x0 = x0_R;
                
                bem.subunits(j).rf_params.left.y0 = y0_L;
                bem.subunits(j).rf_params.right.y0 = y0_R;
            end
            
        end
        
    end
end

function G = gabor(x,y,rf)  
    % Compute two-dimensional vertically oriented Gabor
    % Usage: G = gabor(x,y,rf);
    % x and y: linear arrays; for 2D RFs use meshgrid.
    % rf: RF structure; use rf=bem.default_rf_params() for template
    
    
    G = spatial_gaussian(x,y,rf) .* cos(2*pi*rf.f*x + rf.phi);
end

function G = spatial_gaussian(x,y,rf)
    
    % Compute two-dimensional Gabor
    % x and y: linear arrays, use meshgrid
    % x0 and y0: horizontal and vertical centres, respectively.
    % sx and sy: horizontal and vertical SDs of Gaussian
    G = exp(-((x-rf.x0).^2 ./ (2*rf.sx^2) + ((y-rf.y0).^2)./(2*rf.sy^2)));
end

function G = temporal_gaussian(t,tc,rf)
    tau = rf.tau;

    G = exp(-(t-tc).^2 ./ (2*tau^2));
end

function G = gamma_cosine(t,tc,rf)
    alpha = rf.alpha;
    omega = rf.omega;
    tau = rf.tau;
    t_phi = rf.t_phi;

        G = (1/(gamma(alpha)*tau^alpha) * (t-tc).^(alpha-1) ...
    .* exp(-(t-tc)/tau) .* cos(2*pi*omega*(t-tc) + t_phi))';

    G = real(G); 
    G(t<tc) = 0; 

end

