function tMax = get_maxvar_time(NimStruct,runDur,triggerOnce,windowSize)
    % get_maxvar_time
    % Gets the time with maximum variance across disparities.
    % 
    % Parameters
    % -----------
    % NimCell : NimCell struct
    % 
    % runDur : int, optional
    % Stimulus frame duration (10, 30, 750)
    %
    % triggerOnce : bool, optional
    % Whether to trigger on every frame or just the first frame
    %
    % Returns
    % ---------
    % t : int, time point at which variance across disparity is maximal


    if nargin < 2
        runDur = 30;
    end

    if nargin < 3
        triggerOnce = 1;
    end
    
    if nargin < 4
        windowSize = 15;
    end

    dxs = NimStruct.dxs; %disparity
    corrs = NimStruct.correlation; % binocular correlation 
    durs = NimStruct.duration; % stimulus duration
    dt = NimStruct.times(2)-NimStruct.times(1); % monitor refresh period
    
    runDxs = get_disparities(NimStruct);
    
    if isfield(NimStruct,'rEst'); % this is the case for NimModels
        rObs = NimStruct.rEst;
    else % and this is the case for NimCells
        [rObs,~] = histc(NimStruct.spiketimes,(0:(length(NimStruct.stim)-1))*dt);
    end



    ns = arrayfun(@(k)sum((dxs==k) .* (durs == runDur)),runDxs);
    runDxs = runDxs(ns > max(ns)*0.5);
    
    allResps = zeros(length(runDxs),windowSize);    
    
    for j = 1:length(runDxs)
                
        currentFrames = find((runDxs(j) == dxs) .* (corrs==1).* (runDur == durs) );
            
        if triggerOnce
            currentFrames = remove_duplicates(currentFrames,runDur/10);
        end
        
        rs = get_all_timeseries(rObs,currentFrames);
        
        allResps(j,:) = mean(rs);
        
    end
    
    tVars= var(allResps,[],1);
    tMax = find(tVars==max(tVars));

end

function rs = get_all_timeseries(rObs,currentFrames,windowSize)

    if nargin < 3
        windowSize = 15;
    end

    respCell = arrayfun(@(k)get_timeseries(rObs,k,windowSize),currentFrames,...
        'uniformoutput',0);
    
    rs = cell2mat(respCell)';


end

function y = get_timeseries(rObs,k,windowSize)


    idx = k:(k+windowSize-1);
    idx2 = idx(idx <= length(rObs));
    
    y = zeros(1,windowSize);
    y(1:length(idx2)) = rObs(idx2);
    y =y';
end

function runDxs = get_disparities(NimStruct)
    runDxs = round(unique(NimStruct.dxs),3);   
    runDxs = runDxs(abs(runDxs) < 1e3);
end