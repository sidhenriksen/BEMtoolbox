function [RMS,binocularKernel] = compute_rms(NimModel,binocularKernel)
    % COMPUTE_RMS   computes the RMS contrast for all frames in NimModel.stim.
    % RMS = compute_rms(NimModel,<binocularKernel>)    
    % Computes and returns the RMS for all the frames in NimModel given
    % some model structure (uses the first subunit to estimate RF)
    %
    % Parameters
    % -----------
    % NimModel : NimModel structure
    % binocularKernel (optional) : a binocular kernel to compute the RMS
    % contrast. Default is to use Gaussian estimated from RFs.
        
    if (nargin >= 2)
        
        emptyKernel = isempty(binocularKernel);
        
    else
        
        emptyKernel = false;
        
    end
    
    if (nargin < 2) || emptyKernel
        
        Lx = get_mean_filters(NimModel.NimFit,'L');

        Rx = get_mean_filters(NimModel.NimFit,'R');

        leftGaussian = get_rf_gaussian(Lx);

        rightGaussian = get_rf_gaussian(Rx);

        binocularKernel = [ascolumn(leftGaussian);ascolumn(rightGaussian)];
        
    end
    
    binocularKernel = ascolumn(binocularKernel(:))';

    nFrames = size(NimModel.stim,1);
    
    G = repmat(binocularKernel,[nFrames,1]);
    
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

function Mx = get_filters(NimFit,k,whichEye)

    dims = NimFit.stim_params.dims;    
    dims(2) = dims(2)/2;
    
    
    
    Bvec = NimFit.subunits(k).filtK;

    if strcmp(whichEye,'L')

        Mvec = Bvec(1:length(Bvec)/2);

    elseif strcmp(whichEye,'R')

        Mvec = Bvec((length(Bvec)/2 +1):length(Bvec));

    end

    M = reshape(Mvec',dims(1:2));

    % This is the alternative way of doing it; just takes a cross-section
    [~,tIdx] = max(var(M,[],2));        
    Mx = M(tIdx,:);

    %Mx = std(M,[],1);


end

function Mx = get_mean_filters(NimFit,whichEye)

    dims = NimFit.stim_params.dims;    
    dims(2) = dims(2)/2;
    
    allMx = zeros(dims(2),length(NimFit.subunits));
    for k = 1:length(NimFit.subunits);
        Bvec = NimFit.subunits(k).filtK;

        if strcmp(whichEye,'L')

            Mvec = Bvec(1:length(Bvec)/2);

        elseif strcmp(whichEye,'R')

            Mvec = Bvec((length(Bvec)/2 +1):length(Bvec));

        end

        M = reshape(Mvec',dims(1:2));

        % This is the alternative way of doing it; just takes a cross-section
        %[~,tIdx] = max(var(M,[],2));        
        %Mx = M(tIdx,:);

        allMx(:,k) = std(M,[],1);
    end
    Mx = mean(allMx,2);
end
