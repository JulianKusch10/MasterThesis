% function running_file(taskid)
taskid = 1;
% ncpu = 1;
% maxNumCompThreads(ncpu);

set(0,'defaulttextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

WorkLoc = '.';

mkdir(sprintf('%s/Data',WorkLoc))
mkdir(sprintf('%s/Data/Run_%i',WorkLoc,taskid))

%% Imaginary time
imagtime = 1;
loadHTherm = 1;

%Obtain simulation parameters
[Params] = parameters();

%Set up spatial grids and transforms
[Transf] = setup_space(Params);

%Initialize wavefunction and potentials
[psi,V,VDk] = Initialize(Params,Transf);
%psi = gpuArray(psi);

if loadHTherm == 1
    if Params.T > 0
        load('HThermData.mat',"HThermFit")
        disp('HTherm Loaded')
    elseif Params.T == 0
        kCutoff = 0;
        HThermFit = calculate_Htherm(Params,kCutoff);
        disp('Zero Temperature')
    end
else
    kCutoff = 0;
    HThermFit = calculate_Htherm(Params,kCutoff);
    if Params.T > 0
        save('HThermData.mat',"HThermFit")
    end
     disp('HTherm Calculated')
end

Observ.EVec = []; Observ.NormVec = []; Observ.PCVec = []; Observ.tVecPlot = []; Observ.mucVec = []; 
Observ.res_idx = 1;
t_idx = 1;

[psi] = ssfm_imag(psi,Params,Transf,VDk,V,HThermFit,taskid,t_idx,Observ,WorkLoc);
fprintf('Temperature %i nK calculated',Params.T)