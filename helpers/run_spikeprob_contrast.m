function [binnedContrasts,binnedRates] = run_spikeprob_contrast(NimCell,NimModel,trialMasks,optStruct)
    % Computes the 
    %
    % Usage: [igfs,rates] = run_igf(NimCell,NimFit,trialMasks,optStruct)
    %
    % NimCell : NimCell struct
    % 
    % NimFit : The fitted NIM object (will also take NimModel)
    % 
    % trialMasks : cell array containing logical arrays denoting the
    % triggers for each condition
    %
    % optStruct : optional structure containing fields of whatever 
    % default params to override    
        
        
    nBins = 11; % number of bins
    
    windowSize = 3;
    
    lag = 4;
    
    if nargin > 4
        unpackStruct(optStruct);
    end
        
    
    rObs = get_rate(NimCell);
    %rObs = NimModel.rEst;
    
    rmsContrast = compute_rms(NimModel);
    
    %% Group according to the trialMask
    rmsContrasts = cellfun(@(trialMask) group_data(rmsContrast,trialMask,windowSize,lag),trialMasks,'uniformoutput',false);
    
    rates = cellfun(@(trialMask) group_data(rObs,trialMask,windowSize,lag),trialMasks,'uniformoutput',false);
        
    %bins = prctile(G,linspace(0,1,nBins));
        
    %% Finally bin the databased on igf values
    [binnedContrasts,binnedRates] = cellfun(@(igf,rate)run_binning(igf,rate,nBins),rmsContrasts,rates,'uniformoutput',false);
    
        
end

function [binnedIgf,binnedRate] = run_binning(igf,rate,nBins)

    edges = linspace(0,100,nBins+1);
    edges = edges(1:end-1);
    bins = prctile(igf,edges);
    
    [~,~,idx] = histcounts(igf,bins);    
    
    binnedRate = arrayfun(@(k)(mean(rate(idx==k))),0:(nBins-1));
    binnedIgf = arrayfun(@(k)(mean(igf(idx==k))),0:(nBins-1));

end


function igf = group_data(G,trialMask,windowSize,lag)

    trialIdx = find(trialMask);
    
    igf = zeros(1,length(trialIdx)*windowSize);
    
    for k = 1:length(trialIdx);        
                
        start = trialIdx(k)+lag;
        
        if start > length(trialMask);
            continue
        end
        
        stop = min(start+windowSize-1,length(G));
        
        igfStart = (k-1)*windowSize + 1;
        
        igfStop = min(igfStart+windowSize-1,length(igf));
        
        igf(igfStart:igfStop) = G(start:stop);
        
    end
        

end