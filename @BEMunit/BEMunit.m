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
        tk=struct('tau',[],'omega',[],'alpha',[],'t_phi',[],'t0',[]);
        

        silent = 0;
        
        outputNL = @(x)(x);
        
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
        
        norm_constants=struct();
        
    end
    
    properties (Hidden)
        memory_threshold=5e9; % in bytes; currently 5GB.
        seed_list=[]; % this is for bootstrap_mode
    end
    
    % =======================================================
    % ===== These methods are defined in separate files =====
    methods
        
        bem = add_subunit(bem,varargin);
        bem = modify_subunit(bem,k_sub,varargin);
        bem = update(bem);
        
        [bem,big_bem] = load_bootstrap(bem,generator);
        save_bootstrap(bem,generator);
        delete_bootstrap(bem,generator);
        
        C = simulate_spatial(bem,generator,n_frames,bootstrap_mode,run_parallel,seed);       
        C = simulate_spatiotemporal(bem,generator,n_frames,duration,bootstrap_mode,seed);
        
        
    end
    % =========================================================
    
    
    methods
        
        % Constructor
        function bem = BEMunit(varargin)
            
            % Get the full path of the BEMtoolbox
            BEMunit_dir = mfilename('fullpath');
            stop_idx = strfind(BEMunit_dir,'@BEMunit')-2;
            bem.toolbox_dir=BEMunit_dir(1:stop_idx);
            
            % This is where we save the data
            bem.bootstrap_dir=[bem.toolbox_dir,'/.data/'];   
            
            
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
                    
                    t2 = 0:bem.dt:0.5;
                    
                    rf = bem.subunits(1).rf_params;
                    switch bem.temporal_kernel
                        case 'gaussian'
                            tk=temporal_gaussian(t2,rf.left);                        

                        case 'gamma-cosine'
                            tk=gamma_cosine(t2,rf.left);
                
                    end
                    plot(t2,tk,'k -','linewidth',3);
                    
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
        
        function [cc_dx,ccf] = get_ccf(bem,k_sub)
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
            ccf=conv(xL(end:-1:1),xR);
            cc_idx=round(length(xL)/2+1):round(length(xL)*1.5);
            cc_dx = (cc_idx-median(cc_idx))*bem.deg_per_pixel;
            ccf = ccf(cc_idx);
        end
        
        function bem = compute_normalization_constant(bem,generator,N,bootstrap_mode)
             
             if nargin < 3
                % number of frames to estimate your normalisation constant
                % from
                N = 2500; 
             else
                 if isempty(N);
                     N = 2500;
                 end
             end
             
             if nargin < 4
                 bootstrap_mode = 0;                 
             end

             generator.correlation = 1;
             generator.dx = bem.dx/bem.deg_per_pixel;
             
                          
             if (N == 0) && bootstrap_mode
                 if isfield(bem.subunits(1),'V_L');
                     if isempty(bem.subunits(1).V_L);
                         bem = load_bootstrap(bem,generator);
                     end
                 else
                     bem = load_bootstrap(bem,generator);
                 end
                                  
                 if ~bem.check_bootstrap(generator);
                     N = 2500;
                 end
             end
             
             
             idstring = id2string(generator.get_identifier());
             C = simulate_spatial(bem,generator,N,bootstrap_mode);
             bem.norm_constants.(idstring) = mean(C);             
         end         
         
         function bool = check_normalization_constant(bem,generator);
            generator.correlation = 1;
            generator.dx = bem.dx/bem.deg_per_pixel; 
            idstring = id2string(generator.get_identifier());
            bool = isfield(bem.norm_constants,idstring);
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
            rf_params_both.f=3.125; % cycles per degree
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
       
       function tk = default_tk(~);
           tk.alpha=2.5;
           tk.omega=8.3;
           tk.tau=0.015;
           tk.t_phi=-pi/5;
           tk.t0=0.03;
           
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

            x0_L = bem.x0 - bem.dx/2;
            x0_R = bem.x0 + bem.dx/2;

            y0_L = bem.y0 - bem.dy/2;
            y0_R = bem.y0 + bem.dy/2;
            
            for j = 1:length(bem.subunits);
                bem.subunits(j).rf_params.left.x0 = x0_L;
                bem.subunits(j).rf_params.right.x0 = x0_R;
                
                bem.subunits(j).rf_params.left.y0 = y0_L;
                bem.subunits(j).rf_params.right.y0 = y0_R;
            end
            
        end
            
    end

    methods (Hidden)
        function id = get_identifier(bem,ks,varargin)
            % This method will return a unique identifier
            % based on relevant RF properties for bootstrap resampling.            
            % If 'load_bootstrap' is toggled, this identifier will be used 
            % to find spatial responses that correspond to the cell.
            % Usage: id = bem.get_identifier(<ks>)
            % ks: Optional argument. If not provided, the id will be a
            % unique identifier for the entire unit. If provided, the 
            % id will only be based on the subsets specified
            
            % This is done by adding together the RMS for the left RF,
            % right RF, and the ccf. These should be positionally
            % invariant metrics. We also add the minimum value of
            % each before taking the RMS, so as to make the metric            
            % sensitive to phase information (this is because
            % rms(x)==rms(-x)).
            
            bem_props = {'dim','dx','dy','Nx','Ny','deg_per_pixel'};
            sub_props = {'sx','sy','f','rf_type','phi'};
            
            K = 6;
            
            for k = 1:length(varargin);
                if strcmp(varargin{k},'K');
                    K = varargin{k+1};
                end
            end
            
            all_props = [];
            for prop = 1:length(bem_props);
                current_prop = num2str(bem.(bem_props{prop}));
                all_props = [all_props,current_prop];
            end
            
            if nargin < 2
                ks = 1:bem.n_subunits;
            end
            
            for j = ks
                
                for prop = 1:length(sub_props)                    
                    left_prop = num2str(bem.subunits(j).rf_params.left.(sub_props{prop}));
                    right_prop = num2str(bem.subunits(j).rf_params.right.(sub_props{prop}));
                                        
                    all_props = [all_props,left_prop,right_prop];
                end
            end

            id = stringhash(all_props,'K',K);
            
        end
    

        function bool = check_bootstrap(bem,generator)
            % Method to check whether the current model/stimulus
            % combination has already been bootstrapped
            % Returns a 1 if yes, a 0 if no. 
            % Usage: bool = bem.check_bootstrap(generator);

            bem_id = bem.get_identifier();
            stim_id = generator.get_identifier();

            csvfile = [bem.bootstrap_dir,num2str(bem_id),'_',num2str(stim_id),'.csv'];

            bool=exist(csvfile,'file');

          end
         
    end
    
    methods (Hidden,Static);
        [] = bootstrap_mode_unit_tests();
    end
end
