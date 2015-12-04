classdef RDS < StimulusGenerator
    
    % Default properties
    properties
        dotsize=4; % dot width in pixels
        density=0.25; % dot density; proportion of stimulus coverage assuming no overlap
        dx=0; % horizontal disparity in pixels
        dy=0; % vertical disparity in pixels
        correlation=1; %-1 is anti-correlated; 0 is uncorrelated; 1 is correlated;
        
        Nx=200; % stimulus width in pixels
        Ny=200; % stimulus height in pixels
        
        % 0 = all black dots; 1 = all white dots; 0.5 = equal mix
        polarity=0.5; 

        
        % set to 1 if we want correlation=0 to correspond to half-matched as per Doi et al. (2011).
        toggle_match=0; 
        
    end
    
    properties (Hidden)
        % None
    end
    
    methods
        
        % Constructor
        function rds = RDS(varargin)
            % Constructs an RDS object which can be used to generate
            % binocular stimuli (left and right images).             
            % Usage: rds = RDS(varargin);
            % rds = RDS() will initialise an RDS object with default
            % stimulus parameters
            % 
            % ===optional arguments===
            % dotsize: dot width in pixels
            % density: dot density (fraction of occupied space assuming no
            % dot occlusion; note: dots are allowed to occlude)
            % correlation: binocular correlation of stimulus (-1 to 1)
            % dx: horizontal disparity of stimulus in pixels
            % dy: vertical disparity of stimulus in pixels            
            % Nx: stimulus width in pixels
            % Ny: stimulus height in pixels
            % polarity: polarity of dots (0=all black, 0.5=mixed, 1=all
            % white dots).
            % toggle_match: if 1, then changing correlation will change the
            % *dot match level* as per Doi et al. (2011). Default is 0.
            %
            % Example 1: Generates an anticorrelated RDS at 30% dot
            % density with 2 pixel negative (crossed) disparity.
            % rds = RDS('density',0.3,'correlation',-1,'dx',-2)
            % [L,R] = rds.generate()
            %
            % Example 2: Generates a half-matched RDS with 10px dots
            % rds = RDS('toggle_match',1,'correlation',0,'dotsize',10);
            % [L,R] = rds.generate();
            %
            % Example 3: Properties can also be changed after generator
            % creation. Generates all-white correlated stimuli with 3 pixel
            % negative disparity.
            % rds = RDS();
            % rds.dx = -3;
            % rds.polarity=1;
            % [L,R] = rds.generate();
            
            for j = 1:2:length(varargin);
                assert(isprop(bem,varargin{j}),'Property not in generator structure');
                rds.(varargin{j}) = varargin{j+1};
            end
            
            
        end
        
        function [L,R] = generate(rds)
            % Method to generate a stimulus example.
            % Usage: [L,R] = rds.generate();

            % Parameters for the mother image
            m_Nx = rds.Nx+abs(rds.dx); % x are columns
            m_Ny = rds.Ny+abs(rds.dy); % y are rows
             % dotsize addition needs to be somewhere else

            % The mother images, initiated to gray (0)
            L_mImage = zeros(m_Ny,m_Nx);
            R_mImage = zeros(m_Ny,m_Nx);
            
            dotmatch = (rds.correlation+1)*0.5;
            
            % Number of dots to paint
            n_dots = round(m_Nx*m_Ny * rds.density / (rds.dotsize*rds.dotsize));
            
            if rds.toggle_match
                n_correlated = round(n_dots*dotmatch);
                n_anticorrelated = round(n_dots*(1-dotmatch));
                n_uncorrelated = 0;
            else
                if rds.correlation < 0
                    n_correlated = 0;
                    n_anticorrelated = round(n_dots*abs(rds.correlation));
                else
                    n_correlated = round(n_dots*abs(rds.correlation));
                    n_anticorrelated = 0;
                end
                n_uncorrelated = n_dots - n_correlated - n_anticorrelated;
            end
            
            
            % Add extra centers for uncorrelated because we need to paint those
            % twice (one in each eye).
                        
            % Generate centers; 
            L_centers = [randi(m_Ny-(rds.dotsize-1),n_dots,1),...
                randi(m_Nx-(rds.dotsize-1),n_dots,1)];
            
            R_centers = L_centers; 
            R_centers(n_dots-n_uncorrelated+1:n_dots,:) = ...
                [randi(m_Ny-(rds.dotsize-1),n_uncorrelated,1),...
                randi(m_Nx-(rds.dotsize-1),n_uncorrelated,1)];

            fill_L = [zeros(1,round(n_dots*(1-rds.polarity)))-1, ones(1,round(n_dots*rds.polarity))];
            fill_R = fill_L;
            
            ac_indices = randperm(n_dots-n_uncorrelated,n_anticorrelated); % which ones to anti-correlate
            fill_R(ac_indices) = fill_R(ac_indices)*-1; % invert contrast for ac stimuli
            
            paint_sequence = randperm(n_dots); % randomise the paint sequence
            
            for dot = paint_sequence
                L_mImage(L_centers(dot,1):L_centers(dot,1)+rds.dotsize-1, L_centers(dot,2):L_centers(dot,2)+rds.dotsize-1) = fill_L(dot);
                R_mImage(R_centers(dot,1):R_centers(dot,1)+rds.dotsize-1, R_centers(dot,2):R_centers(dot,2)+rds.dotsize-1) = fill_R(dot);
            end


            if rds.dx <= 0 && rds.dy <= 0
                L = L_mImage(1:rds.Ny,1:rds.Nx);
                R = R_mImage((1:rds.Ny)-rds.dy, (1:rds.Nx)-rds.dx);

            elseif rds.dx <= 0 && rds.dy > 0
                L = L_mImage((1:rds.Ny)+rds.dy,1:rds.Nx);
                R = R_mImage(1:rds.Ny, (1:rds.Nx)-rds.dx);

            elseif rds.dx > 0 && rds.dy > 0
                L = L_mImage(1:rds.Ny,(1:rds.Nx)+rds.dx);
                R = R_mImage((1:rds.Ny)+rds.dy, 1:rds.Nx);

            elseif rds.dx > 0 && rds.dy <= 0
                L = L_mImage(1:rds.Ny,(1:rds.Nx)+rds.dx);
                R = R_mImage((1:rds.Ny)-rds.dy, 1:rds.Nx);
            end
        end
    end
    
end