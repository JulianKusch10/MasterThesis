function [psi] = ssfm_imag(psi,Params,Transf,VDk,V,HThermFit,taskid,t_idx,Observ,WorkLoc)

set(0,'defaulttextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

dt=-1j*abs(Params.dt);

KEop= 0.5*(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);
Observ.residual = 1e-2; Observ.res = 1; 

HT = HTherm(HThermFit,Params,psi);
muchem = chemicalpotential(psi,Params,Transf,VDk,V,HT);

figure(1)
runningplot(psi,Params,Transf,Observ)
drawnow

%
Ftherm = EthermInterp(Params,HThermFit);  %construct the thermal energy interpolant
%

AdaptIdx = 0;
psi = real(psi);

while t_idx < Params.cut_off
    %kin
    psi = fftn(psi);
    psi = psi.*exp(-0.5*1i*dt*KEop);
    psi = ifftn(psi);
    
    %DDI
    Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
    psi = sqrt(Params.N)*psi/sqrt(Norm);
    frho = fftn(abs(psi).^2);
    Phi = ifftn(frho.*VDk);

    %Real-space
    HT = HTherm(HThermFit,Params,psi);
    psi = psi.*exp(-1i*dt*(V + Params.gs*abs(psi).^2 + Params.gammaQF*abs(psi).^3 + Params.gdd*Phi + HT - muchem));
    % 
    % Norm = sum(abs(psi1(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
    % psi1 = sqrt(Params.N)*psi1/sqrt(Norm);
    % 
    % frho = fftn(abs(psi1).^2);
    % Phi = real(ifftn(frho.*VDk));
    % psi = psi.*exp(-1i*dt*(V + Params.gs*abs(psi1).^2 + Params.gammaQF*abs(psi1).^3 + Params.gdd*Phi - muchem));

    %kin
    psi = fftn(psi);
    psi = psi.*exp(-0.5*1i*dt*KEop);
    psi = ifftn(psi);
   
    % %Renorm
    Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
    psi = sqrt(Params.N)*psi/sqrt(Norm);

    muchem = chemicalpotential(psi,Params,Transf,VDk,V,HT);

    %Plotting loop
    if mod(t_idx,1000) == 0

        %Change in Energy
        E = energytotal(psi,Params,Transf,VDk,V,Ftherm);
        E = E/Norm;
        Observ.EVec = [Observ.EVec E];

        %Chemical potential
        Observ.mucVec = [Observ.mucVec muchem];

        %Normalized residuals
        res = norm_resid(psi,Params,Transf,VDk,V,muchem,HT);
        Observ.residual = [Observ.residual res];
        
        Observ.res_idx = Observ.res_idx + 1;
        figure(1)
        runningplot(psi,Params,Transf,Observ)
        drawnow

        save(sprintf('%s/Data/Run_%i/psi_gs.mat',WorkLoc,taskid),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');

        %Adaptive time step -- Careful, this can quickly get out of control
        relres = abs(Observ.residual(Observ.res_idx)-Observ.residual(Observ.res_idx-1))/Observ.residual(Observ.res_idx);
        if relres < 1e-4
            if AdaptIdx > 2 && abs(dt) > Params.mindt
                dt = dt / 2;
                fprintf('Time step changed to '); disp(dt);
                AdaptIdx = 0;
            elseif AdaptIdx > 3 && abs(dt) < Params.mindt
                break
            else
                AdaptIdx = AdaptIdx + 1;
            end
        else
            AdaptIdx = 0;
        end
        % psi = fftn(psi);
    end
    if any(isnan(psi(:)))
        disp('Idiot.')
        break
    end
    t_idx=t_idx+1;
end

%Change in Energy
E = energytotal(psi,Params,Transf,VDk,V,Ftherm);
E = E/Norm;
Observ.EVec = [Observ.EVec E];

%Chemical potential
Observ.mucVec = [Observ.mucVec muchem];

%Normalized residuals
res = norm_resid(psi,Params,Transf,VDk,V,muchem,HT);
Observ.residual = [Observ.residual res];

Observ.res_idx = Observ.res_idx + 1;

save(sprintf('%s/Data/Run_%i/psi_gs.mat',WorkLoc,taskid),'psi','muchem','Observ','t_idx','Transf','Params','VDk','V');
