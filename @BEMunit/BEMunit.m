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
        n_subunits=0;
        
        init=true;
        
        bootstrap_dir;
        bootstrap_loaded=false;
        toolbox_dir;
        
    end
    
    properties (Hidden)
        memory_threshold=5e9; % in bytes; currently 5GB.
    end
    
    methods
        % These methods are defined in separate files
        bem = add_subunit(bem,varargin);
        bem = modify_subunit(bem,k_sub,varargin);
        bem = update(bem);
        bem = load_bootstrap(bem);
        
        C = simulate_spatial(bem,generator,n_frames);
        C = simulate_spatiotemporal(bem,generator,n_frames,duration,bootstrap_mode);
        
        
    end
    
    methods
        
        % Constructor
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
                
            end

            % Get the full path of the BEMtoolbox
            BEMunit_dir = mfilename('fullpath');
            stop_idx = strfind(BEMunit_dir,'@BEMunit')-2;
            bem.toolbox_dir=BEMunit_dir(1:stop_idx);
            
            % This is where we save the data
            bem.bootstrap_dir=[bem.toolbox_dir,'/.data/'];            
            
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
                        [cc_idx,ccf] = get_ccf(bem,j);
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
        
        function [cc_idx,ccf] = get_ccf(bem,k_sub)
            % Get cross-correlation function of the kth subunit.
            % Usage: ccf = bem.get_ccf(k);
            % k: index of subunit
            % Returns the cross-sectional cross-correlation function of 
            % the left and right receptive field.
            
            L = bem.subunits(k_sub).L;
            R = bem.subunits(k_sub).R;
            [~,max_idx] = max(var(L,[],1));

            % cross section of left and right RF
            xL = L(max_idx,:);
            xR = R(max_idx,:);

            % cross-correlation function
            ccf=conv(xL,xR(end:-1:1));
            cc_idx=round(length(xL)/2+1):round(length(xL)*1.5);
            ccf = ccf(cc_idx);                        
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
                bem.x0 = (rand-0.5)*0.1*x_max;
            end
            if isempty(bem.y0);
                bem.y0 = (rand-0.5)*0.1*y_max;
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
        
                
        function id = get_identifier(bem,ks)
            % This method will return a unique identifier
            % based on relevant RF properties for bootstrap resampling.            
            % If 'load_bootstrap' is toggled, this identifier will be used 
            % to find spatial responses that correspond to the cell.
            % Usage: id = bem.get_identifier(<ks>)
            % ks: Optional argument. If not provided, the id will be a
            % unique identifier for the entire unit. If provided, the 
            % 
            
            % This is done by adding together the RMS for the left RF,
            % right RF, and the ccf. These should be positionally
            % invariant metrics. We also add the minimum value of
            % each before taking the RMS, so as to make the metric            
            % sensitive to phase information (this is because
            % rms(x)==rms(-x)).
            
            id=0;
            
            if nargin ~= 2
                ks = 1:bem.n_subunits;
            end
               
            
            for j = ks
                L = bem.subunits(j).L;
                R = bem.subunits(j).R;
                ccf = get_ccf(bem,j); % cross-correlation function
                
                id = id+rms(L(:)+min(L(:))) + rms(R(:)+min(R(:))) + rms(ccf+min(ccf(:)));
            end
            
        end
        
        function save_bootstrap(bem)
            
    
        end
        
        
        
    end
end
