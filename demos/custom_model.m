

bem = BEMunit('init',0); % this will not initialise to a standard energy model

% Creates a subunit with phase pi/2 in the left eye [L_phi], 0 in the right eye [R_phi],
% and whose orientation preference is pi/8 on both eyes [L_theta and R_theta].
bem = bem.add_subunit('L_phi',pi/2,'R_phi',0,'L_theta',pi/8,'R_theta',pi/8);

% Creates a subunit with a half-wave rectification followed by a squaring
% (default is squaring)
bem = bem.add_subunit();
bem.subunits(2).outputNL = @(x)(((x>0).*x).^2);

bem.dx = 0.15;

bem = bem.update();

bem.display_model()

% Stimulus parameters
rds = RDS();
N = 1e3;
dxs = -10:15;

tc = zeros(1,length(dxs));
for j = 1:length(dxs);
    rds.dx = dxs(j);
    tc(j) = mean(bem.simulate_spatial(rds,N));
end

figure();
plot(dxs*bem.deg_per_pixel,tc./max(tc),'k o -','linewidth',3,'markersize',8);
xlabel('Disparity (deg)','fontsize',16);
ylabel('Normalized response','fontsize',16);
set(gcf,'color','white');
set(gca,'fontsize',14);
