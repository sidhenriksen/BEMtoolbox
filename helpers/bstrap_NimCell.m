function bstrapNimCell = bstrap_NimCell(Xcell,rObs,p)
    % this function will resample the frames of NimCell
    % with replacement. The rationale here is that the
    % trials are not independent, and so resampling might
    % be able to get decent performance while reducing the
    % duration of training.
    %
    % === Parameters ===
    % NimCell : 
    % p : number of samples, expressed as a proportion of total N of trials
    %
    % === Returns ===
    % bstrap_NimCell
    
    N = length(rObs);
    Nresample = round(N*p);
    
    idx = randi(N,Nresample,1);
    
    simpleFields = {'dxs','times','correlation','duration'}

end