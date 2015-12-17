% SetupBEMtoolbox();

%%% Compute BEM responses to a dynamic RDS

bem = BEMunit();
bem.Nx = 125; bem.Ny=125;
bem.deg_per_pixel=0.04;
bem.dx=0.12;
% default is sx=0.1; rescale to sx=0.15; also changes freq of gabor
bem = bem.rescale(1.5);

bem.temporal_kernel = 'gaussian';
bem = bem.update(); % update bem structure after making changes

% Display model:
bem.display_model()

rds = RDS();
rds.Nx = bem.Nx; rds.Ny = bem.Ny;
rds.density=0.8; % 80% density stimuli
rds.dotsize=3;
rds.dx=bem.dx/bem.deg_per_pixel;

duration = 2; % 2s duration
n_frames = 40; % 40 frames; pattern refresh rate is 20Hz

rds.correlation = 0.5;
rds.display_stimulus(); title('Correlation = 0.5');
C1 = bem.simulate_spatiotemporal(rds,n_frames,duration);


rds.correlation = -0.5;
rds.display_stimulus(); title('Correlation=-0.5')
C2 = bem.simulate_spatiotemporal(rds,n_frames,duration);

t = (0:(length(C1)-1))*bem.dt;
figure(); hold on;
plot(t,C1,'k -','linewidth',2);
plot(t,C2,'r -','linewidth',2)
legend('Correlation=0.5','Correlation=-0.5');
xlabel('Time (s)','fontsize',16);
ylabel('Response');