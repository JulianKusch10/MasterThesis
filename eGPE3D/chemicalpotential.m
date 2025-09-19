function muchem = chemicalpotential(psi,Params,Transf,VDk,V)

%Parameters
KEop= 0.5*(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);

% DDIs
frho=fftn(abs(psi).^2);
Phi=real(ifftn(frho.*VDk));

Eddi = (Params.gdd*Phi.*abs(psi).^2);
Eddi = sum(Eddi(:))*Transf.dx*Transf.dy*Transf.dz;

%Kinetic energy
% Ekin = KEop.*abs(fftn(psi)*normfac).^2;
% Ekin = sum(Ekin(:))*Transf.dkx*Transf.dky*Transf.dkz/(2*pi)^3;

Ekin = conj(psi).*ifftn(KEop.*fftn(psi));
Ekin = sum(Ekin(:))*Transf.dx*Transf.dy*Transf.dz;

%Potential energy
Epot = V.*abs(psi).^2;
Epot = sum(Epot(:))*Transf.dx*Transf.dy*Transf.dz;

%Contact interactions
Eint = Params.gs*abs(psi).^4;
Eint = sum(Eint(:))*Transf.dx*Transf.dy*Transf.dz;

%Quantum fluctuations
Eqf = Params.gammaQF*abs(psi).^5;
Eqf = sum(Eqf(:))*Transf.dx*Transf.dy*Transf.dz;

%Total energy
muchem = Ekin + Epot + Eint + Eddi + Eqf; %
muchem = real(muchem) / Params.N;