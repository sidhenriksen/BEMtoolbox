function NimCell = make_NimCell(AllStims,AllSpikes,cellNo)
    % Function to take allStims and allSpikes structure and cast it to a
    % stim and spiketimes structure used by the NIM toolbox
    
    NUM_FRAMES_PER_TRIAL=300;
    
    % This gives a stimIndex and spikeIndex pair which matches trials in
    % allStims with the spike data
    [stimIndices,spikeIndices] = get_trial_indices(AllStims,AllSpikes,cellNo);
    
    
    % Get the data for the current cell
    currentData = AllSpikes(cellNo);
    datafields = fields(currentData);
    for k = 1:length(datafields);
        df = datafields{k};
        if length(currentData.(df)) == length(spikeIndices);
            currentData.(df) = currentData.(df)(spikeIndices);
        end
    end
        
    spikeData = {currentData.spikes{spikeIndices}};    
    stimData = AllStims(stimIndices);
    
    
    
    dxData = currentData.dx;
    corrData = currentData.corrs;

    
    Ns = zeros(1,length(stimData));
    for j = 1:length(stimData)
        Ns(j) = size(stimData(j).L,2);
    end
    nmin = min(Ns);
    
    for j = 1:length(stimData);
        stimData(j).L = stimData(j).L(1:NUM_FRAMES_PER_TRIAL,:);
        stimData(j).R = stimData(j).R(1:NUM_FRAMES_PER_TRIAL,:);
        stimData(j).contrast = stimData(j).contrast(1:NUM_FRAMES_PER_TRIAL,:);
        
        if size(stimData(j).L,2) > nmin;
            stimData(j).L = stimData(j).L(:,1:nmin);
            stimData(j).R = stimData(j).R(:,1:nmin);
        end
        
        dxData{j} = round(convert_dx(dxData{j}),3);
        corrData{j} = convert_dx(corrData{j});
    end
    
    
    L = cat(1,stimData.L);
    R = cat(1,stimData.R);
    contrast = cat(1,stimData.contrast);
    
    stim = [L,R];
    
    
    dxs = [];    
    corr = [];    
    dur = [];
    
    
    keepIndices = [];
    startIndices = 1;
    
    
    spiketimes = [];
    times = [];
            
    relativeSpiketimes = [];
    trials = [];
    
    

    trialStart = zeros(1,length(spikeData));
    trialEnd = zeros(1,length(spikeData));
    % So this puts all the spiketimes into a somewhat fake absolute time.
    % Also want to discard trials with less than, say 20 frames
    k = 0;
    N = zeros(1,length(spikeData));
    
    for trial = 1:length(spikeData)
        
        nFrames = size(stimData(trial).L,1);
        
        N(trial) = nFrames;
        
        currentDx = dxData{trial};
        currentCorr = corrData{trial};
        
        
        stopIndices = startIndices+(nFrames-1);
        
        % stimStart and stimStop track which stimuli we want to include in
        % the final segment
        if nFrames >= 20;
            k = k+1; % Only increment counter if we're including this trial
            if k==1;
                trialStart(k) = 0;
            else
                trialStart(k) = trialEnd(k-1) + 0.01;
            end

            trialEnd(k) = trialStart(k) + (nFrames-1)*0.01;

            currentSpikes = double(spikeData{k}')/10000;
            currentSpiketimes = currentSpikes(currentSpikes>0) + trialStart(k);

            relativeSpiketimes = [relativeSpiketimes,currentSpikes(currentSpikes>0)];
            spiketimes = [spiketimes,currentSpiketimes];

            addTimes = (0:(nFrames-1))*0.01;
            times = [times,addTimes];
            dxs = [dxs,currentDx(1:nFrames)];
            corr = [corr,currentCorr(1:nFrames)];
            dur = [dur,zeros(1,nFrames)+currentData.dur(trial)];

            nSpikes = length(currentSpiketimes);
            trials = [trials,zeros(1,nSpikes)+k];

            keepIndices = [keepIndices,startIndices:stopIndices];


            
        end
        startIndices = stopIndices+1;
        
         
        
    end
    
    stim = stim(keepIndices,:);
    contrast = contrast(keepIndices,:);
    
    NimCell.stim = stim;
    NimCell.spiketimes = spiketimes;
    NimCell.dxs = dxs;
    NimCell.times = times;
    NimCell.trials = trials;
    NimCell.correlation = corr;
    NimCell.duration = dur;
    NimCell.contrast = contrast;
    
    
   
end

function y = convert_dx(x)

    K = 300;
    
    if iscolumn(x)
        x = x';
    end
    
    if length(x) < K;
        rep = K/length(x);
        
        y = repmat(x,[rep,1]);
        y = y(:)';
        y = y(1:K);
    else
        y = x;
    end

end