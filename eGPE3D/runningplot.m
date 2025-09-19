function RPlot = runningplot(psi,Params,Transf,Observ)
    set(0,'defaulttextInterpreter','latex')
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    
    format long
%     percentcomplete = linspace(0,1,Params.cut_off/200);
   
    %Plotting    

    subplot(2,3,1)
    n = abs(psi).^2;
    nxz = squeeze(sum(n*Transf.dy,2));
    nyz = squeeze(sum(n*Transf.dx,1));
    nxy = squeeze(sum(n*Transf.dz,3));
    
    plotxz = pcolor(Transf.x,Transf.z,nxz');
    set(plotxz, 'EdgeColor', 'none');
    xlabel('$x$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');
       
    subplot(2,3,2)
    plotyz = pcolor(Transf.y,Transf.z,nyz');
    set(plotyz, 'EdgeColor', 'none');
    xlabel('$y$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');

    subplot(2,3,3)
    plotxy = pcolor(Transf.x,Transf.y,nxy');
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
