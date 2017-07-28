function varargout = plot_tc(tcMatrix,fig,lineStyle,x,lw)
    % Takes as input a 3 x N matrix, with 3 correlations,
    % and N disparities.
    % Usage: 
    % plot_tc(tcMatrix,handle,lineStyle,dx,lineWidth);
    % 
    

    if nargin < 2 || isempty(fig)
        fig = figure();
        hold on; 
    end
    
    if nargin < 3
        lineStyle = '-';
    end
    
    if nargin < 4
        x = 1:size(tcMatrix,2);
    end
    
    if nargin < 5
        lw = 3;
    end
    
    if strfind(class(fig),'axis')
        ax=fig;
    else
        figure(fig);
        ax=gca;
    end
    
    
    plot(ax,x,tcMatrix(1,:),'k','linewidth',lw,'linestyle',lineStyle)
    plot(ax,x,tcMatrix(2,:),'color',[0.7,0.7,0.7],'linestyle',lineStyle,'linewidth',lw);
    plot(ax,x,tcMatrix(3,:),'r','linewidth',lw,'linestyle',lineStyle);
    
    if nargout
        varargout = {fig};
    end
end