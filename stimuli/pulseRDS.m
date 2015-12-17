classdef pulseRDS < StimulusGenerator
    
    
    properties    
        pwm=10;    
        status=0; % init to off
        count=0; % count is 0
        generator;
    end
    
    % constructor
    methods
        function rds = pulseRDS(generator)
            rds.generator=generator;
        end
        
        function [L,R] = generate(rds)
            rds.count = rds.count+1;
            if rds.count >= rds.pwm;
                rds.status = ~rds.status;
                rds.count = 0;
            end
            
            [L,R] = rds.generator.generate();
            L=L*rds.status;
            R=R*rds.status;
            
        end
    end
    
    
    
end