function fig3()
    
    red = [0.8,0.1,0.1];
    blue = [0.1,0.1,0.8];

    constants = figure_constants();
        
    [allCells] = get_plotter_data();
    
    m = arrayfun(@(K)mean(strip_uc(K.twopassSpikeCount)),allCells);
    v = arrayfun(@(K)mean(strip_uc(K.totalVariance)),allCells);
    ev = arrayfun(@(K)mean(strip_uc(K.totalVariance-K.internalVariance)),allCells);
    
    figure();
    subplot(1,3,1); hold on;    
    plot(m,v,'o','markersize',constants.markersize,'markerfacecolor',red,'color',red);
    xlabel('Mean spike count','fontsize',constants.fontsize)
    ylabel('Total variance','fontsize',constants.fontsize)
    
    x = linspace(0,max(m)*1.1,500);
    y1 = x + x.^2/3;
    y2 = x + x.^2;
    
    plot(x,y1,'--','linewidth',constants.linewidth,'color',blue)
    plot(x,y2,'--','linewidth',constants.linewidth,'color',blue)
    
    ylim([0,5])
    
    %title('Total variance','fontsize',constants.fontsize)
    
    
    subplot(1,3,2); hold on;
    
    plot(m,ev,'o','markersize',constants.markersize,'markerfacecolor',red,'color',red);
    xlabel('Mean spike count','fontsize',constants.fontsize)
    ylabel('External variance','fontsize',constants.fontsize)
    
    x = linspace(0,max(m)*1.1,500);
    y1 = x.^2/3;
    y2 = x.^2;
    
    plot(x,y1,'--','linewidth',constants.linewidth,'color',blue)
    plot(x,y2,'--','linewidth',constants.linewidth,'color',blue)
    
    ylim([0,3])
    
    %title('External variance','fontsize',constants.fontsize)
        
    subplot(1,3,3); hold on;
    
    plot(m,ev./m,'o','markersize',constants.markersize,'markerfacecolor',red,'color',red);
    xlabel('Mean spike count','fontsize',constants.fontsize)
    ylabel('External Fano Factor','fontsize',constants.fontsize)
    
    x = linspace(0,max(m)*1.1,500);
    y1 = x;
    y2 = x/3;
    
    plot(x,y1,'--','linewidth',constants.linewidth,'color',blue)
    plot(x,y2,'--','linewidth',constants.linewidth,'color',blue)
    
    ylim([0,4])
    
    %title('Fano Factor','fontsize',constants.fontsize);
    
    set(gcf,'color','white')

    for k = 1:3;
        
        subplot(1,3,k);
        
        set(gca,'fontsize',constants.fontsize)
        
    end 
end