function [lowIdx,highIdx] = split_trials_by_contrast(NimModel,p,runDur)

    if nargin < 2
        p = 1/3;
    end
    
    if nargin < 3
        runDur = 30;
    end
    
    idx = NimModel.duration == runDur;
    
    Lx = get_filters(NimModel.NimFit,'L');
    
    Rx = get_filters(NimModel.NimFit,'R');
    
    leftGaussian = get_rf_gaussian(Lx);
    
    rightGaussian = get_rf_gaussian(Rx);
    
    binocularGaussian = [leftGaussian,rightGaussian];
    
    RMS = compute_rms(NimModel,binocularGaussian);
    
    lowIdx = RMS < quantile(RMS(idx),p);
    
    highIdx = RMS > quantile(RMS(idx),1-p);
     
end

function RMS = compute_rms(NimModel,binocularGaussian);

    N = size(NimModel.stim,1);
    
    G = repmat(binocularGaussian,[N,1]);
    
    RMS = rms(G.*NimModel.stim,2);
    
end

function G = get_rf_gaussian(M)
                    
    Mrect = abs(M);
            
    Mrect = Mrect-median(Mrect);
        
    x = 1:length(M);
    
    params = fit_gaussian(x,Mrect);
    
    % set offset to 0 and amplitude to 1
    params(1) = 0; params(2) = 1;
        
    
              
    G = get_gaussian(params,x);
        
end

function Mx = get_filters(NimFit,whichEye)

    dims = NimFit.stim_params.dims;    
    dims(2) = dims(2)/2;
    
    Bvec = NimFit.subunits(1).filtK;
    
    if strcmp(whichEye,'L')
        
        Mvec = Bvec(1:length(Bvec)/2);
        
    elseif strcmp(whichEye,'R')
        
        Mvec = Bvec((length(Bvec)/2 +1):length(Bvec));
        
    end
        
    M = reshape(Mvec',dims(1:2));
    
    % This is the alternative way of doing it; just takes a cross-section
    %[~,tIdx] = max(var(M,[],2));        
    %Mx = M(tIdx,:);
        
    Mx = std(M,[],1);
        
end
