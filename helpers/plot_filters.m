
function plot_filters(NimModel)

    useSameLims = 0;

    if isfield(NimModel,'NimFit')
        
        NimFit = NimModel.NimFit;
        
    else
        
        NimFit = NimModel;
        
    end
    
    allFilts = [NimFit.subunits.filtK];
    
    cMin = min(allFilts(:));
    cMax = max(allFilts(:));
    
    dims = NimFit.stim_params.dims(1:2);
    
    K = length(NimFit.subunits);
    [nRows,nCols] = get_subplot_dims(K);
    
    figure();
    for k = 1:K
        subplot(nRows,nCols,k);
        
        if NimFit.subunits(k).weight == 1
            suType = 'excitatory';
        else
            suType = 'suppressive';
        end
        
        
        if ~useSameLims
            cMin = min(allFilts(:,k));
            cMax = max(allFilts(:,k));
            cAbs = max([abs(cMin),abs(cMax)]);
            cMin = -cAbs;
            cMax = cAbs;
        end

        B = reshape(NimFit.subunits(k).filtK,dims);
        imagesc(B);
        set(gca,'clim',[cMin,cMax]);
        title(sprintf('%s',suType),'fontsize',12);

    end
    

end

function [nRows,nCols] = get_subplot_dims(K)
    % returns subplot dimensions for K neurons
    
    while isprime(K) && (K>3)
        K = K+1;
    end
    
    
    costs = zeros(K,K);
    for row = 1:K
        for col = 1:K
            
            if row*col == K
                costs(row,col) = abs(row-col);
            else
                costs(row,col) = inf;
            end
        end
    end
    
    [~,idx] = min(costs(:));
    
    [nRows,nCols] = ind2sub(size(costs),idx);
            
    

end