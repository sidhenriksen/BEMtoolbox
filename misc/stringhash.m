function hash=stringhash(S,varargin)
    % Quick and shitty hash function with probably very terrible
    % collision probability. If bootstrap_mode fails it is almost 
    % guaranteed to be because of this.
    % Usage: hash = stringhash(S)
    % S: A string
    % Will return a number between 0 and 2^32.
        
    assert(ischar(S),'Error: Input must be a string');
    
    M = 2^32;   
    
    ks = linspace(1,2,length(S));
    S2 = round(S.^ks *100);
    hash = mod(sum(S2),M);
end