function run_variance_regression()

    [allCells,allModels] = get_plotter_data();

    
    R1s = zeros(1,length(allCells));
    R2s = zeros(1,length(allCells));
    R3s = zeros(1,length(allCells));

    for k = 1:length(allCells)

        externalVarianceCell = strip_uc(allCells(k).totalVariance - allCells(k).internalVariance);

        externalVarianceModel = strip_uc(allModels(k).totalVariance - allModels(k).internalVariance);
        
        spikeCountCell = strip_uc(allCells(k).twopassSpikeCount);
        
        spikeCountModel = strip_uc(allModels(k).twopassSpikeCount);        
        
        meanCountDifference = ascolumn(spikeCountCell-spikeCountModel);
        meanSquaredCountDifference = ascolumn(spikeCountCell.^2-spikeCountModel.^2);
        
        y = ascolumn(externalVarianceCell)-ascolumn(externalVarianceModel);
        
        X1 = meanCountDifference;
        
        X2 = [X1,meanSquaredCountDifference];
        
        noOffset = 0;
        
        r1 = regression2(y,X1,noOffset);
        
        noOffset = 1;
        
        r2 = regression2(y,X1,noOffset);
        
        r3 = special_regression(X1,y);
        
        R1s(k) = r1;
        R2s(k) = r2;
        R3s(k) = r3;

    end
    
    bins = linspace(0,1,21);
    figure();
    histogram(R1s.^2,bins,'facecolor',[0.8,0.1,0.1],'facealpha',1);
    xlabel('R^2','fontsize',16);
    ylabel('Frequency','fontsize',16);
    xlim([-0.05,1.05]);
    ylim([0,9]);
    set(gca,'fontsize',13);
    
    title(sprintf('Mean=%.2f, median=%.2f, SD=%.2f',mean(R1s.^2),median(R1s.^2),std(R1s.^2)),'fontsize',16);
    
    
    

end

function r = special_regression(x,y)

    cost = @(a,b)sum((a-b).^2);
    
    fitfunc = @(x,y,params)(cost(x*params,y));
    
    init = 1;
    
    myAns = fminsearch(@(params)fitfunc(x,y,params),init);
    
    r = corrcoef(ascolumn(x*myAns),ascolumn(y));
    r = r(2);

end


    
