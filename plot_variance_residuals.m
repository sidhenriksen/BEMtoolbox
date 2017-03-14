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
    externalVarianceCell = (allCells(k).totalVariance - allCells(k).internalVariance);

    externalVarianceModel = (allModels(k).totalVariance - allModels(k).internalVariance);

    spikeCountCell = strip_uc(allCells(k).twopassSpikeCount);

    spikeCountModel = strip_uc(allModels(k).twopassSpikeCount);        

    meanCountDifference = ascolumn(spikeCountCell-spikeCountModel);
    %meanSquaredCountDifference = ascolumn(spikeCountCell.^2-spikeCountModel.^2);

    
    
    
    y = ascolumn(strip_uc(externalVarianceCell - externalVarianceModel));

    X1 = meanCountDifference;
    
    noOffset = 1;
    [r,m,b] = regression2(y,X1,noOffset); 
    
    xrange = range(X1);
    t = linspace(min(X1)-xrange*0.1,max(X1)+xrange*0.1,51);
    g = b + m*t;

    externalVarianceRegModel = restore_uc(X1*m) + externalVarianceModel;

    allValues = [externalVarianceCell(:);externalVarianceModel(:);externalVarianceRegModel(:)];
    
    lims = [min(allValues(:)),max(allValues(:))];
    
    ax1 = subplot(2,2,1);
    cla; hold on;
    plot_tc(externalVarianceCell,ax1);
    title('External variance (cell)','fontsize',14);
    ylim(lims);
    
    ax2 = subplot(2,2,2);
    cla; hold on;
    plot_tc(externalVarianceModel,ax2);
    title('External variance (model)','fontsize',14);
    ylim(lims);
    
    ax3 = subplot(2,2,3);
    cla; hold on;
    plot_tc(externalVarianceRegModel,ax3);
    title('External variance (regression model)','fontsize',14);
    ylim(lims);
    
    subplot(2,2,4);
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


function y = restore_uc(x)


    N = (length(x)-1)/2;
    
    uc = x(end);
    x = x(1:end-1);
        
    y(1,:) = x(1:2:end);
    y(2,:) = uc;
    y(3,:) = x(2:2:end);

end
