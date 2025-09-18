    set(0,'defaulttextInterpreter','latex')
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    format long

	runIdx = 2;
	WorkLoc = 'E:/ComplexLangevin/NoLHY';
            %'/gpfs/bwfor/work/ws/hd_tf297-BdG3D';
            %'E:/ComplexLangevin';
            %'.';
	load(sprintf('%s/Data/Run_%i/psi_gs.mat',WorkLoc,runIdx),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');
    
    x = Transf.x*Params.l0*1e6; 
    y = Transf.y*Params.l0*1e6; 
    z = Transf.z*Params.l0*1e6; 
%     percentcomplete = linspace(0,1,Params.cut_off/200);
    
    dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
    %Plotting    
    subplot(2,3,1)
    n = abs(psi).^2;
    nxz = squeeze(sum(n*dy,2));
    nyz = squeeze(sum(n*dx,1));
    nxy = squeeze(sum(n*dz,3));
    
    plotxz = pcolor(x,z,nxz');
    set(plotxz, 'EdgeColor', 'none');
    xlabel('$x$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');
       
    subplot(2,3,2)
    plotyz = pcolor(y,z,nyz');
    set(plotyz, 'EdgeColor', 'none');
    xlabel('$y$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');

    subplot(2,3,3)
    plotxy = pcolor(x,y,nxy');
    set(plotxy, 'EdgeColor', 'none');
    xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]');

    subplot(2,3,4)
    plot(-log10(Observ.residual),'-b')
    ylabel('$-\mathrm{log}_{10}(r)$'); xlabel('steps');

    subplot(2,3,5)
    plot(Observ.EVec,'-b')
    ylabel('$E$'); xlabel('steps');

    subplot(2,3,6)
    plot(Observ.mucVec,'-b')
    ylabel('$\mu$'); xlabel('steps');
%     xlim([0,1]); ylim([0,8]);
%     xlim([0,1]); ylim([0,8]);

	Ecomp = energy_components(psi,Params,Transf,VDk,V);
