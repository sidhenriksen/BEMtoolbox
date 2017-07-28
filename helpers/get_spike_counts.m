function spikeCounts = get_spike_counts(rEst,idx,lag,step)
    % GET_SPIKE_COUNTS  returns the spike count for each frame in idx
    % spikeCounts = get_spike_counts(rEst,idx,<lag>,<step>)
    % 
    % Parameters
    % -----------
    % rEst : array of spikes
    % idx : index array (logical or by normal index)
    % lag (optional) : lag after idx to sum. Default is 4 frames.
    % step (optional) : step across which to sum. Default is 3 frames.
    %
    % Returns
    % --------
    % spikeCounts : array of spike counts, with the same dimensions as idx
    % (or smaller if idx+lag+step exceeds bounds).
    
    if nargin < 4
        step = 3;
    end
    
    if nargin < 3
        lag = 4;
    end
    
    
    if islogical(idx)
        idx = find(idx);
    end
    
    if isempty(idx)
        spikeCounts = []; return;                
    end
    
    baseIdx = repmat(ascolumn(idx),[1,step]);
    
    addFrames = repmat( (0:step-1) +lag,[length(idx),1]);
        
    newIdx = baseIdx + addFrames;
    
    % this deletes any frames that go on for longer than the length
    % of rEst.
    popIndex = newIdx(end-30:end,end) > length(rEst);
    popIndex = logical([zeros(size(newIdx,1)-31,1);popIndex]);
    newIdx(popIndex,:) = [];
    
    
    if size(newIdx,2) > 1
        spikeCounts = sum(rEst(newIdx),2);
    else
        spikeCounts = rEst(newIdx);
    end

end