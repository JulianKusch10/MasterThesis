% function running_file(taskid)
taskid = 1;
% ncpu = 1;
% maxNumCompThreads(ncpu);

WorkLoc = '.';
            %'/gpfs/bwfor/work/ws/hd_tf297-BdG3D';
            %'.';

mkdir(sprintf('%s/Data',WorkLoc))
mkdir(sprintf('%s/Data/Run_%i',WorkLoc,taskid))

%% Switches
imagtime_switch = 1;
quench_switch = 0;

%% Imaginary time

if imagtime_switch == 1
    %Obtain simulation parameters
    [Params] = parameters(taskid);
    
    %Set up spatial grids and transforms
    [Transf] = setup_space(Params);
    
    %Initialize wavefunction and potentials
    [psi,V,VDk] = Initialize(Params,Transf);
    %psi = gpuArray(psi);

    Observ.EVec = []; Observ.NormVec = []; Observ.PCVec = []; Observ.tVecPlot = []; Observ.mucVec = []; 
    Observ.res_idx = 1;
    t_idx = 1;
    [psi] = ssfm_imag(psi,Params,Transf,VDk,V,taskid,t_idx,Observ,WorkLoc);
    disp('Imaginary time completed')
else
    load(sprintf('%s/Data/Run_%i/psi_gs.mat',WorkLoc,taskid),'psi','Params','Transf','VDk','V','muchem')
    disp('Initial state loaded from file')
end 

%% === Dynamics === %%
if quench_switch == 1
    
    mkdir(sprintf('%s/Data/Run_%i/TimeEvolution',WorkLoc,taskid))
    Params.dt = 0.0002;
    Quench = QuenchSettings(Params, taskid);
    save(sprintf('%s/Data/Run_%i/QuenchSettings',WorkLoc,taskid),'Quench','Params');
    psi0 = psi;
    
    if Quench.addnoise == 1
        rng('shuffle')
        [noise] = thermalnoise3D(Transf,Params,Quench.Temp,0);
        psi = psi0 + noise;
        Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz; % normalisation
        psi = sqrt(Params.N)*psi/sqrt(Norm);
    end
        
    % --- Initialize
    Observ.EVec = []; Observ.NormVec = []; Observ.PCVec = []; Observ.tVecPlot = []; Observ.asVec = [];
    Observ.res_idx = 1; 
    [psi] = ssfm_real(psi,psi0,Params,Transf,VDk,V,taskid,Observ,WorkLoc,Quench);

end