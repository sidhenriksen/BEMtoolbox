function [lowIdx,highIdx] = split_trials_by_contrast(NimModel,binocularKernel,optStruct)
    % SPLIT_TRIALS_BY_CONTRAST  splits the NimModel frames into low and high contrast.
    % [lowIdx,highIdx] =  split_trials_by_contrast(NimModel,<p>,<runDur>,<binocularKernel>,<pTrials>)    
    %
    % Parameters
    % -----------
    % NimModel : a NimModel structure
    % binocularKernel (optional) : a custom binocular kernel. Gaussian fit
    % is default.
    % 
    % optStruct (optional) : a structure with fields corresponding to any
    % of the less commonly changed parameters.
    %   p : determines centile split. p=1/3 is default.
    %   runDur : frame duration to use. 30 ms is default.
    %   pTrials : proportion of trials to use. 0.8 is default.
    %   cv : boolean specifying to run on test set (cv=0) or CV set (cv=1).
    %   Default is 0.
    %   seed : this is random 
    %        
    % Returns
    % -------
    % lowIdx : logical array where true entries correspond to low contrast
    % frames
    % highIdx : logical array where true entries correspond to high contrast
    % frames
    
    p = 1/3;
        
    %pTrials = 0.8;
    
    cv = 0;
    
    seed = 'shuffle';

    %runDur = 30;

    if nargin > 3

        unpackStruct(optStruct);
    end
    
    if strcmp(seed,'shuffle') && cv
        
        warning('Cross-validation is enabled, but seed is set to random.')
        
    end
 
    rng(seed);

    

    %idx = NimModel.duration == runDur;
    
    %flipidx = flipbit(idx,1-pTrials);
    
     if cv
         NimModel.stim = NimModel.stimCv;
         idx = NimModel.indexCv;

%         idx = idx & ~flipidx;
%         
     else
         NimModel.stim = NimModel.stimTrain;
         idx = NimModel.indexTrain;

%         idx = idx & flipidx;
     end

    % NimModel.stim = NimModel.stim(idx,:);
    
    if nargin < 2
        RMS = compute_rms(NimModel);
    elseif isempty(binocularKernel)
        RMS = compute_rms(NimModel);
    else
        RMS = compute_rms(NimModel,binocularKernel);
    end
    
    % RMS2 is to account for the fact that we have used only a subsample
    % of the trials . In order to quickly return a logical index, we need
    % to create a dummy array with NaNs filled in where we didn't compute
    % RMS
    
    RMS2 = zeros(size(idx));
    RMS2(idx) = RMS;
    RMS2(~idx) = NaN;
    
    
    lowIdx = idx & (RMS2 < quantile(RMS,p));
    
    highIdx = idx & (RMS2 > quantile(RMS,1-p));
     
end

function y = flipbit(x,p);

    N = round(length(x)*p);
    
    flipIdx = randperm(length(x),N);
    
    y = x; 
    
    y(flipIdx) = ~y(flipIdx);

end


