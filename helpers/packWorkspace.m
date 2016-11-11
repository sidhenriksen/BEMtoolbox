% This function packs all workspace variables into a struct.
%
% If some arguments were supplied, only those variables whose names match
% an argument will be packed.
%
% Ghaith Tarawneh (ghaith.tarawneh@ncl.ac.uk)
%
function struct = packWorkspace(varargin)

workspaceVars = evalin('caller', 'who');

for i=1:length(workspaceVars)
    
    w = workspaceVars{i};
    
     if nargin > 0 % if vararin is supplied
        
        if sum(strcmp(varargin, w)) == 0
            
            continue; % current field not in varargin, ignore
            
        end
        
    end
    
    thisvar = evalin('caller', w);
    
    struct.(workspaceVars{i}) = thisvar;
end

