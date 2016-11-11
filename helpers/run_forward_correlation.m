function [runDxs,tc,vars,allResps] = run_forward_correlation(NimStruct,runDur,windowSize,lag,triggerOnce,trialMask)
    % Runs forward correlation on NimStruct to construct disparity tuning
    % curves for different binocular correlations.
    %
    % Usage: [runDxs,tc,vars,allResps] = ...
    % run_forward_correlation(NimCell,<runDur>,<windowSize>,<lag>,<triggerOnce>,<trialMask>)
    %
    % NimStruct : NimCell or NimModel struct
    % 
    % runDur (optional) : Frame duration to trigger frames on (default 30 ms)
    % windowSize (optional) : size of the forward window to run forward
    % correlation; at 10 ms resolution, windowSize=15 would mean 150 ms
    % window).
    %
    % lag (optional) : An optional delay (e.g., 2 would delay the beginning
    % of the window to 20 ms after the frame onset)
    %
    % triggerOnce : Whether to trigger on only the first frame presentation,
    % e.g. the first frame of a sequence of 3x10 ms frames.
    %
    % trialMask : logical array specifying which trials to include 
    %
    
    
    if nargin < 2
        runDur = 30;
    end
    
    if nargin < 3
        windowSize = 15;
    end
    
    if nargin < 4
        lag = 0;
    end
    
    if nargin < 5
        triggerOnce = 1;
    end
    
    if nargin < 6
        trialMask = ones(size(NimStruct.dxs));
    end
    
    trialMask = ascolumn(trialMask)';
    
    dxs = round(NimStruct.dxs,3); %disparity
    corrs = NimStruct.correlation; % binocular correlation 
    durs = NimStruct.duration; % stimulus duration
    dt = NimStruct.times(2)-NimStruct.times(1); % monitor refresh period
    
    runDxs = get_disparities(NimStruct);
    
    runCorrs = [-1,0,1];
    
    if isfield(NimStruct,'rEst'); % this is the case for NimModels
        rObs = NimStruct.rEst;
    else % and this is the case for NimCells
        [rObs,~] = histc(NimStruct.spiketimes,(0:(length(NimStruct.stim)-1))*dt);
    end
        
    % initialise arrays
    tc = zeros(length(runCorrs),length(runDxs));
    vars = zeros(length(runCorrs),length(runDxs));
    allResps = cell(length(runCorrs),length(runDxs));
    
    for j = 1:length(runDxs);
        for i = 1:length(runCorrs);
            
            if runCorrs(i) ~= 0;
                
                currentFrames = find((runDxs(j) == dxs) .* (corrs==runCorrs(i)).* (runDur == durs) .* ...
                    trialMask);
            else
                currentFrames = find((corrs==runCorrs(i)).* (runDur == durs) .* trialMask);
            end
            
            if triggerOnce
                currentFrames = remove_duplicates(currentFrames,runDur/10);
            end

            resps = zeros(1,length(currentFrames));
            
            for k = 1:length(currentFrames);
                start = currentFrames(k)+lag;
                stop = min([start+windowSize-1,length(dxs)]);

                R = zeros(windowSize,1);

                
                R(1:length(start:stop)) = rObs(start:stop);

                resps(k) = sum(rObs(start:stop));

            end

            tc(i,j) = mean(resps);
            vars(i,j) = var(resps);
            
            allResps{i,j} = resps;
        end
    end
    
    Ns = cellfun(@length,allResps(1,:));
    include = Ns./mean(Ns) > 0.25; % greater than 25% of the mean N    
    runDxs = runDxs(include);
    tc = tc(:,include);
    vars = vars(:,include);
    allResps = allResps(:,include);
    
    if sum(~include) > 0;
        fprintf('Excluded %i disparities.\n',sum(~include))
    end
end

function runDxs = get_disparities(NimStruct)
    runDxs = round(unique(NimStruct.dxs),3);   
    runDxs = runDxs(abs(runDxs) < 1e3);
end