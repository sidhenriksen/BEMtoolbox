classdef StimulusGenerator < handle
   
    properties
        
    end
    
    properties (Hidden)
        
    end
    
    methods
        % Constructor
        function generator = StimulusGenerator()
            
        end
        
        
        function id = get_identifier(generator)

            props = properties(generator);
            prop_id = [];

            for i = 1:length(props);
                prop_value = generator.(props{i});
                prop_id = [prop_id,props{i},num2str(prop_value)];
            end
            
            big_id = [class(generator),prop_id];
            id=stringhash(big_id);
            
        end
        
        function unit_tests(generator,N)
            % Function will run some unit tests using whichever properties
            % are available for the current generator.
            % Usage: generator.unit_tests(...)
            % Optional input:
            % N: Number of repeats of the stimulus (default is 25).
                        
            
            if nargin < 2;
                N = 25;
            end
            
            % These are here just to avoid any errors due to calling properties
            % that do not exist.
            if isprop(generator,'correlation');
                correlation = generator.correlation;
            else
                correlation = 1;
            end
            
            if isprop(generator,'dx');
                dx = generator.dx;
            else
                dx = 0;
            end
            
            ds = zeros(1,N);
            cs = zeros(1,N);
            lums = zeros(2,N);
            dxs = zeros(1,N);
            
            M = zeros(generator.Nx*generator.Ny,1);
            
            for j = 1:N;
                
                % Compute density
                [L,R] = generator.generate();
                M = M+abs(L(:))+abs(R(:));
                % Density mean (sd)
                d = (mean(L(:)~=0) + mean(R(:)~=0))/2;
                ds(j) = d;

                
                
                % Compute binocular correlation
                % set dx to 0 so that we can compute correlation easily
                if isprop(generator,'dx')
                    generator.dx = 0;
                end
                [L0,R0] = generator.generate();                                

                % stick correlation values in array; reset dx value 
                c = corr(L0(:),R0(:));
                cs(j) = c;
                if isprop(generator,'dx')
                    generator.dx = dx;
                end
                
                

                % Compute mean luminance (left and right)
                lums(1,:) = mean(L(:));
                lums(2,:) = mean(R(:));
                
                
                
                % Set disparity to proper value; set correlation to 1
                % Estimate disparity:
                if isprop(generator,'correlation')
                    generator.correlation = 1;
                    
                end
                
                [Lc,Rc] = generator.generate();

                xL = Lc(round(size(Lc,1)/2),:);
                xR = Rc(round(size(Rc,1)/2),:);

                % cross-correlation function
                ccf=conv(xL,xR(end:-1:1));
                [~,peak_idx] = max(ccf);
                dx_est = length(xL)-peak_idx;
                dxs(j) = dx_est;

                if isprop(generator,'correlation')
                    generator.correlation = correlation;
                end
                
            end
            
            % This is a makeshift correction for non-plottable area in a
            % stimulus (for example gray-space where there will never be
            % any dots. These should not be counted in the density estimate.
            density_normalisation = sum(M(:)~=0)/length(M);
            ds = ds/density_normalisation;
            
            %results_table = [mean(ds),std(ds),median(ds)
            
            fprintf('Metric       \t Mean  \t\t SD  \t\t Median  \t Real\n')
            if isprop(generator,'density')
                fprintf('Density      \t %.3f  \t %.3f  \t %.3f  \t %.3f\n',mean(ds),std(ds),median(ds),generator.density)
            else
                fprintf('Density      \t %.3f  \t %.3f  \t %.3f  \t N/A\n',mean(ds),std(ds),median(ds))
            end
            if isprop(generator,'correlation')
                fprintf('Correlation  \t %.3f  \t %.3f  \t %.3f  \t %.3f\n',mean(cs),std(cs),median(cs),generator.correlation)
            else
                fprintf('Correlation  \t %.3f  \t %.3f  \t %.3f  \t N/A\n',mean(cs),std(cs),median(cs))
            end
            
            fprintf('Luminance L  \t %.3f  \t %.3f  \t %.3f  \t N/A\n',mean(lums(1,:)),std(lums(1,:)),median(lums(1,:)))
            fprintf('Luminance R  \t %.3f  \t %.3f  \t %.3f  \t N/A\n',mean(lums(2,:)),std(lums(2,:)),median(lums(2,:)))
            
            if isprop(generator,'dx')
                fprintf('Estimated dx \t %.3f  \t %.3f  \t %.3f  \t %.3f\n',mean(dxs),std(dxs),median(dxs),generator.dx)
            else
                fprintf('Estimated dx \t %.3f  \t %.3f  \t %.3f  \t N/A\n',mean(dxs),std(dxs),median(dxs))
            end
            
            
            
        end
        
        function display_stimulus(generator);
            [L,R] = generator.generate();
            cmin = min([L(:);R(:)]);
            cmax = max([L(:);R(:)]);
            
            figure();
            colormap('gray')
            subplot(1,2,1);
            imagesc(L);
            set(gca,'clim',[cmin,cmax],'xtick',[],'ytick',[]);
            
            subplot(1,2,2);
            imagesc(R);
            set(gca,'clim',[cmin,cmax],'xtick',[],'ytick',[]);
            
        end
        
        function newGenerator = copy(generator)
           
            newGenerator = feval(class(generator));
            
            p = properties(generator);
            for i = 1:length(p)
                newGenerator.(p{i}) = generator.(p{i});
            end
            
        end
        
    end
    
end