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
    
    
    
    
        
%     
%     
%     K = 4;
%     
%     
%     for j = 1:length(varargin);
%         if strcmp(varargin{j},'K');
%             K = varargin{j+1};
%         end
%         
%         if strcmp(varargin{j},'M');
%             M = varargin{j+1};
%         end
%     end
%     
%     N = ceil(length(S)/K);
%     
%     mults = (0:(K-1)).*256;
%     mults(1) = 1;
%     
%     hash = 0;
%     for j = 1:N;
%         start = (j-1)*K+1;
%         stop = (j-1)*K + K;
%         stop = min(length(S),stop);
%         
%         current_S = S(start:stop);
%         hash = hash + current_S*mults(1:length(current_S))';
%         
%     end
% 
%     hash = mod(abs(hash),M);   
end