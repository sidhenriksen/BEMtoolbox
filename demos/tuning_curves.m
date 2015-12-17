
% SetupBEMtoolbox(); 

%%% Compute a basic disparity tuning curve

% Set model parameters
bem = BEMunit();
bem.Nx = 100; bem.Ny = 100; %dimensions for model/stimulus
bem.deg_per_pixel = 0.025;
bem.dx = -0.1; % set disparity
bem = bem.update(); % this updates the model with any changes you've made

% Set stimulus parameters
rds = RDS();
rds.Nx = bem.Nx; rds.Ny = bem.Ny;
dxs = (-12:12) + bem.dx/bem.deg_per_pixel;

tcs = zeros(2,length(dxs)); % array to hold the tuning curve

N = 1000; % How many RDSs to average across

for j = 1:length(dxs);
    rds.dx = dxs(j);
    
    rds.correlation = -1;    
    AC = bem.simulate_spatial(rds,N);
    tcs(1,j) = mean(AC);
    
    rds.correlation = 1;
    C = bem.simulate_spatial(rds,N);
    tcs(2,j) = mean(C);
    
end

figure(); hold on;
plot(dxs*bem.deg_per_pixel,tcs(1,:)./max(tcs(:)),'k o --','linewidth',3,'markersize',8)
plot(dxs*bem.deg_per_pixel,tcs(2,:)./max(tcs(:)),'r o --','linewidth',3,'markersize',8)
ylim([0,1]);
xlim([-0.45,0.25])
xlabel('Disparity (deg)','fontsize',16);
ylabel('Response','fontsize',16)
set(gcf,'color','white');
set(gca,'fontsize',14);