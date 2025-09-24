function [psi] = ssfm_real(psi, Params, Transf, VDk, V, taskid, Observ, WorkLoc, Quench, HThermFit)

set(0,'defaulttextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%n0 = abs(psi0).^2;
%dn_ik = fftn(abs(psi).^2-n0);

KEop= 0.5*(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);
t_idx = 1;
dt = Params.dt;
%J_z = zeros(length(Transf.z),Quench.ObservSteps);
%Skz = zeros(length(Transf.z),Quench.ObservSteps);

%
%Normalization
Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
Observ.NormVec = [Norm];

Ftherm = EthermInterp(Params,HThermFit);  %construct the thermal energy interpolant

%Change in Energy
Observ.res_idx = 1;
E = energytotal(psi, Params, Transf, VDk, V, Ftherm);
E = E/Norm;
Observ.EVec = [E];

%Chemical potential
HT = HTherm(HThermFit,Params,psi);
muchem = chemicalpotential(psi, Params, Transf, VDk, V, HT);
Observ.mucVec = [muchem];



%Current density J_z(t)
%[~,~,psi_z] = gradient(psi,Transf.dx,Transf.dy,Transf.dz);
%[~,~,psi_z_T] = gradient(conj(psi),Transf.dx,Transf.dy,Transf.dz);
%J_full = -1j*0.5*(conj(psi).*psi_z - psi.*psi_z_T);
%J_z(:,Observ.res_idx) = J_full(end/2+1,end/2+1,:);

%Density structure factor S(kz,t)
%dn_tk = fftn(abs(psi).^2-n0);
%dndn = fftshift(dn_tk.*dn_ik);
%Skz(:,Observ.res_idx) = ifftshift(dndn(end/2+1,end/2+1,:)); %ifftshift should prep for the 2D fftshift in the analysis (time direction)

%Time
Observ.tVecPlot = 0;

%
%figure(1)
%runningplot_realtime(psi,Params,Transf,Observ,Quench)
%drawnow
% save(sprintf('%s/Data/Run_%i/TimeEvolution/psi_%i.mat',WorkLoc,99,Observ.res_idx),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');

h = waitbar(0,'Please wait...');
while t_idx <= Quench.tSteps
    current_time = t_idx * Quench.dt_ms;
    waitbar(t_idx/Quench.tSteps, h, sprintf('Step %d of %d or %.3f ms of %.0f ms',t_idx ,Quench.tSteps, current_time, Quench.total_time));
    [~, ~, VDk] = Initialize(Params,Transf);
    % Parameters at time t
    tVal = Quench.tVec(t_idx);
    Params.as = Quench.as_vec(t_idx);
    Params.gs = Quench.gs_vec(t_idx);
    Params.gammaQF = Quench.gammaQF_vec(t_idx);
    Params.T = Quench.Temp_vec(t_idx);

    %kin
    psi = fftn(psi);
    psi = psi.*exp(-0.5*1i*dt*KEop);
    psi = ifftn(psi);
    
    %DDI
    frho = fftn(abs(psi).^2);
    Phi = real(ifftn(frho.*VDk));

    %Real-space
    psi = psi.*exp(-1i*dt*(V + Params.gs*abs(psi).^2 + Params.gammaQF*abs(psi).^3 + Params.gdd*Phi - muchem));    

    %kin
    psi = fftn(psi);
    psi = psi.*exp(-0.5*1i*dt*KEop);
    psi = ifftn(psi);
   
    %Chemical potential
    muchem = chemicalpotential(psi, Params, Transf, VDk, V, HT);

    %Plotting loop
    if mod(t_idx,Quench.savestep) == 0
        Observ.res_idx = Observ.res_idx + 1;

        %Current density J_z(t)
        %[~,~,psi_z] = gradient(psi,Transf.dx,Transf.dy,Transf.dz);
        %[~,~,psi_z_T] = gradient(conj(psi),Transf.dx,Transf.dy,Transf.dz);
        %J_full = -1j*0.5*(conj(psi).*psi_z - psi.*psi_z_T);
        %J_z(:,Observ.res_idx) = J_full(end/2+1,end/2+1,:);
%
        %%Density structure factor S(kz,t)
        %dn_tk = fftn(abs(psi).^2-n0);
        %dndn = fftshift(dn_tk.*dn_ik);
        %Skz(:,Observ.res_idx) = ifftshift(dndn(end/2+1,end/2+1,:)); %ifftshift should prep for the 2D fftshift in the analysis (time direction)

        %Normalization
        Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
        Observ.NormVec = [Observ.NormVec Norm];
        
        %Change in Energy
        E = energytotal(psi, Params, Transf, VDk, V, Ftherm);
        E = E/Norm;
        Observ.EVec = [Observ.EVec E];

        %Chemical potential
        Observ.mucVec = [Observ.mucVec muchem];
  	
	    %time
	    Observ.tVecPlot = [Observ.tVecPlot tVal];
        title = sprintf('Step %d of %d or %.3f ms of %.0f ms',t_idx ,Quench.tSteps, current_time, Quench.total_time);
        figure(1)
        runningplot_realtime(psi, Params, Transf, Observ, Quench, title)
        drawnow

        % save(sprintf('%s/Data/Run_%i/TimeEvolution/psi_%i.mat',WorkLoc,99,Observ.res_idx),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');
    end

    if any(isnan(psi(:)))
        disp('Idiot.')
        break
    end
    t_idx=t_idx+1;
end

Observ.res_idx = Observ.res_idx + 1;

%Normalization
Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
Observ.NormVec = [Norm];

%Change in Energy
E = energytotal(psi, Params, Transf, VDk, V, Ftherm);
E = E/Norm;
Observ.EVec = [E];

%Chemical potential
muchem = chemicalpotential(psi, Params, Transf, VDk, V, HT);
Observ.mucVec = [muchem];

%Current density J_z(t)
%[~,~,psi_z] = gradient(psi,Transf.dx,Transf.dy,Transf.dz);
%[~,~,psi_z_T] = gradient(conj(psi),Transf.dx,Transf.dy,Transf.dz);
%J_full = -1j*0.5*(conj(psi).*psi_z - psi.*psi_z_T);
%J_z(:,Observ.res_idx) = J_full(end/2+1,end/2+1,:);
%
%%Density structure factor S(kz,t)
%dn_tk = fftn(abs(psi).^2-n0);
%dndn = fftshift(dn_tk.*dn_ik);
%Skz(:,Observ.res_idx) = ifftshift(dndn(end/2+1,end/2+1,:)); %ifftshift should prep for the 2D fftshift in the analysis (time direction)

%time
Observ.tVecPlot = [Observ.tVecPlot tVal];

% save(sprintf('%s/Data/Run_%i/TimeEvolution/psi_%i.mat',WorkLoc,99,Observ.res_idx),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');
%save(sprintf('%s/Data/Run_%i/TimeEvolution/Jz_%i.mat',WorkLoc,taskid,QuenchIdx),'J_z','Observ','t_idx','Transf','Params');
%save(sprintf('%s/Data/Run_%i/TimeEvolution/Skz_%i.mat',WorkLoc,taskid,QuenchIdx),'Skz','Observ','t_idx','Transf','Params');