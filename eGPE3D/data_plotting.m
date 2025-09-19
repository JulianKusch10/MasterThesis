set(0,'defaulttextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

y = load("Data\Run_1\psi_gs.mat");
n = abs(psi).^2;
nxy = squeeze(sum(n*Transf.dz,3));

clf;
plotxy = pcolor(Transf.x,Transf.y,nxy');
set(plotxy, 'EdgeColor', 'none');
set(gca, 'Color', 'white');
xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]');
