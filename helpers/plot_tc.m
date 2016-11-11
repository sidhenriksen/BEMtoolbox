function varargout = plot_tc(tcMatrix,fig,lineStyle);
    % Takes as input a 3 x N matrix, with 3 correlations,
    % and N disparities.


    if nargin < 2
        fig = figure();
        hold on; 
    end
    
    if nargin < 3
        lineStyle = '-';
    end
    
    x = 1:size(tcMatrix,2);
    
    if strfind(class(fig),'axis')
        ax=fig;
    else
        figure(fig);
        ax=gca;
    end
    
    
    plot(ax,x,tcMatrix(1,:),'k','linewidth',2,'linestyle',lineStyle)
    plot(ax,x,tcMatrix(2,:),'color',[0.7,0.7,0.7],'linestyle',lineStyle);
    plot(ax,x,tcMatrix(3,:),'r','linewidth',2,'linestyle',lineStyle);
    
    if nargout
        varargout = {fig};
    end
end