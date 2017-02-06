
[allCells,allModels] = get_plotter_data();

k = 20;

dx = allCells(k).dx;
dx = [dx,dx(end)+(dx(2)-dx(1))];    

m = allCells(k).twopassSpikeCount;

v = allCells(k).totalVariance;

lw = 3;

figure();

subplot(1,3,1); hold on;
plot(dx,m(1,:),'k -','linewidth',lw);
plot(dx,m(2,:),'-','linewidth',lw,'color',[0.6,0.6,0.6]);
plot(dx,m(3,:),'r -','linewidth',lw);
%title('Mean spike count','fontsize',16);
xlabel('Disparity (deg)','fontsize',16)
ylabel('Mean spike count','fontsize',16)

subplot(1,3,2); hold on;
plot(dx,v(1,:),'k -','linewidth',lw);
plot(dx,v(2,:),'-','linewidth',lw,'color',[0.6,0.6,0.6]);
plot(dx,v(3,:),'r -','linewidth',lw);
xlabel('Disparity (deg)','fontsize',16)
ylabel('Spike count variance','fontsize',16)

ff = v./m;
subplot(1,3,3); hold on;
plot(dx,ff(1,:),'k -','linewidth',lw);
plot(dx,ff(2,:),'-','linewidth',lw,'color',[0.6,0.6,0.6]);
plot(dx,ff(3,:),'r -','linewidth',lw);
xlabel('Disparity (deg)','fontsize',16)
ylabel('Fano Factor','fontsize',16)
ylim([1,2])

for k = 1:3
    subplot(1,3,k);
    set(gca,'fontsize',14,'linewidth',2);
    
end

set(gcf,'color','white')