% This is a little helper function that will "unpack" a structure; it will
% load all the fields as independent variables in the workspace of the
% function that calls it.
%
% This can be useful in scenarios where all default experiment parameters
% are coded within the script and you need an easy way to override these by
% passing a structure containing new values to the function, for example:
%
% function foobar(expt)
% a = 1;
% b = 2;
% c = 3;
% if nargin > 0
%   loadExpt(expt)
% end
% % do some work
% end
%
% In the above, the values of a, b and c can be overridden by passing them
% as fields in the structure expt. If no struct is passed, the default
% values are used (1, 2, 3 in the example).
%
% In addition to expt, the user may supply a list of variable names. If
% this list is supplied, only variables whose names are in the list will be
% unpacked.
%
% Ghaith Tarawneh (ghaith.tarawneh@ncl.ac.uk)
%
function unpackStruct(expt, varargin)

if ~isstruct(expt)
    
    error('Argument is not a struct.')
    
end

fns = fieldnames(expt);

vals = struct2cell(expt);

l = length(fns);

for i=1:l

    fn = fns{i};
    
    fv = vals{i};
    
    if nargin > 1 % if vararin is supplied
        
        if sum(strcmp(varargin, fn)) == 0
            
            continue; % current field not in varargin, ignore
            
        end
        
    end
    
    assignin('caller', fn, fv);
    
end

end