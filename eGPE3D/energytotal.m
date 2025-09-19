function E = energytotal(psi,Params,Transf,VDk,V)

%Parameters

KEop= 0.5*(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);
normfac = Params.Lx*Params.Ly*Params.Lz/numel(psi);

% DDIs
frho = fftn(abs(psi).^2);
Phi = ifftn(frho.*VDk);

Eddi = 0.5*Params.gdd*Phi.*abs(psi).^2;
Eddi = sum(Eddi(:))*Transf.dx*Transf.dy*Transf.dz;

%Kinetic energy
% Ekin = KEop.*abs(fftn(psi)*normfac).^2;
% Ekin = sum(Ekin(:))*Transf.dkx*Transf.dky*Transf.dkz/(2*pi)^3;
% 
Ekin = conj(psi).*ifftn(KEop.*fftn(psi));
Ekin = sum(Ekin(:))*Transf.dx*Transf.dy*Transf.dz;

% Potential energy
Epot = V.*abs(psi).^2;
Epot = sum(Epot(:))*Transf.dx*Transf.dy*Transf.dz;

%Contact interactions
Eint = 0.5*Params.gs*abs(psi).^4;
Eint = sum(Eint(:))*Transf.dx*Transf.dy*Transf.dz;

%Quantum fluctuations
Eqf = 0.4*Params.gammaQF*abs(psi).^5;
Eqf = sum(Eqf(:))*Transf.dx*Transf.dy*Transf.dz;

%
E = Ekin + Epot + Eint + Eddi + Eqf;
E = real(E);