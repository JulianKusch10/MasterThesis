function res = norm_resid(psi,Params,Transf,VDk,V,muchem)

KEop= 0.5*(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);

% DDIs
frho=fftn(abs(psi).^2);
Phi=ifftn(frho.*VDk);

Eddi = Params.gdd*Phi.*psi;

%Kinetic energy
Ekin = ifftn(KEop.*fftn(psi));
Ekin = Ekin;

%Potential energy
Epot = V.*psi;

%Contact interactions
Eint = Params.gs*abs(psi).^2.*psi;

%Quantum fluctuations
Eqf = Params.gammaQF*abs(psi).^3.*psi;

%Total energy
res =  sum(abs(Ekin(:) + Epot(:) + Eint(:) + Eddi(:) + Eqf(:) - muchem*psi(:))*Transf.dx*Transf.dy*Transf.dz)/sum(abs(muchem*psi(:))*Transf.dx*Transf.dy*Transf.dz);