

function make_plain(fig)

    ms = 14;

    ch = get(fig,'children');
    
    ax = ch(end);
    
    axch = get(ax,'children');
    
   
    for k = 1:length(axch);

            set(axch(k),'markersize',ms,'marker','o');

        
    end

    title([]);
    
    set(gca,'fontsize',14);
end
