classdef StimulusGenerator
   
    properties
        
    end
    
    properties (Hidden)
        
    end
    
    methods
        % Constructor
        function generator = StimulusGenerator()
            
        end
        
        
        function id = get_id(generator)
            id = sum(class(generator));
            props = properties(generator);
            
            for i = 1:length(props);
                
                id = id + 1./(1+abs(sum(props{i})));
                
            end
            
        end
        
        %function [L,R] = generate(generator)
            % Generate method for generating an instance of the stimulus
            
        %end
        
    end
    
end