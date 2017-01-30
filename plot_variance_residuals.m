function plot_variance_residuals(allCells,allModels)

    if ~nargin
        [allCells,allModels] = get_plotter_data();
    end
    
    fig=figure();
    set(fig,'windowkeypressfcn',@update_figure);
    
    setappdata(gcf,'allCells',allCells);
    setappdata(gcf,'allModels',allModels);
    setappdata(gcf,'ctr',1);
    
    plot_data(allCells,allModels,1);

end

function plot_data(allCells,allModels,k)
    externalVarianceCell = strip_uc(allCells(k).totalVariance - allCells(k).internalVariance);

    externalVarianceModel = strip_uc(allModels(k).totalVariance - allModels(k).internalVariance);

    spikeCountCell = strip_uc(allCells(k).twopassSpikeCount);

    spikeCountModel = strip_uc(allModels(k).twopassSpikeCount);        

    meanCountDifference = ascolumn(spikeCountCell-spikeCountModel);
    %meanSquaredCountDifference = ascolumn(spikeCountCell.^2-spikeCountModel.^2);

    y = ascolumn(externalVarianceCell)-ascolumn(externalVarianceModel);

    X1 = meanCountDifference;
    
    [r,m,b] = regression2(y,X1);        
    
    xrange = range(X1);
    t = linspace(min(X1)-xrange*0.1,max(X1)+xrange*0.1,51);
    g = b + m*t;

    plot(X1,y,'k o','markersize',7,'markerfacecolor','k'); hold on;
    plot(t,g,'r -','linewidth',2);
    hold off;
    xlabel('Mean spike count difference'); ylabel('External variance residuals');

    cellName = allCells(k).fileName(1:end-4);
    title(sprintf('%s, R^2=%.2f',cellName,r^2),'fontsize',16);
    drawnow;
end

function update_figure(fig,evt)    

    
    allCells = getappdata(fig,'allCells');
    
    allModels = getappdata(fig,'allModels');
    
    k = getappdata(fig,'ctr');
    
    switch evt.Key
        case 'leftarrow'
            if k ~= 1;
                k = k-1;
            end
        case 'rightarrow'
            if k ~= length(allCells)
                k = k+1;
            end
    end
    
    plot_data(allCells,allModels,k);
    setappdata(fig,'ctr',k);
            
    

end
