function VDk = VDcutoff(kr, kz, Rmax, Zmax, Nr)
% makes the dipolar direct interaction matrix, size numel(kr) * numel(kz)
% Rmax and Zmax are the interaction cutoffs. I use 4*(where the density goes to 10^-4 of its peak)
% VDk needs to be multiplied by Cdd
% approach is that of Lu, PRA 82, 023622 (2010)
% blame Danny Baillie, 9 Aug 2011

% accuracy inputs for numerical integration
if(nargin==4)
    Nr = 5e4;
end
Nz = 64;
farRmultiple = 2000;

% analytical transform without cutoff
[KR, KZ]=ndgrid(kr,kz);
Ksq = KR.^2 + KZ.^2;
cossq = KZ.^2./Ksq;
VDk = cossq-1/3;

% analytical cutoff for slice 0<z<Zmax, 0<r<Inf Ronen, PRL 98, 030406 (2007)
sinsq = 1 - cossq;
VDk = VDk + exp(-Zmax*KR).*( sinsq .* cos(Zmax * KZ) - sqrt(sinsq.*cossq).*sin(Zmax * KZ) );

% midpoint grids for the integration over 0<z<Zmax, Rmax<r<farRmultiple*Rmax (i.e. starts at Rmax)
dr=(farRmultiple-1)*Rmax/Nr;
r = ((1:Nr)'-0.5)*dr+Rmax;
dz=Zmax/Nz;
z = ((1:Nz)-0.5)*dz;
[R, Z] = ndgrid(r,z);

Rsq = R.^2 + Z.^2;

% real space interaction to be transformed
igrandbase = (1 - 3*Z.^2./Rsq)./Rsq.^(3/2);

% do the Fourier transform numerically
% prestore to ensure each besselj is only calculated once
% cell is faster than (:,:,krn) slicing
Nkr = numel(kr);
besselr = cell(Nkr,1);
 for krn = 1:Nkr
    besselr{krn} = repmat(r.*besselj(0,kr(krn)*r),1,Nz);
 end
 for kzn = 1:numel(kz) % what goes wrong when kzn = 33?
    igrandbasez = repmat(cos(kz(kzn)*z),Nr,1) .* igrandbase;
    for krn = 1:Nkr
        igrand = igrandbasez.*besselr{krn};
        VDk(krn,kzn) = VDk(krn,kzn) - sum(igrand(:))*dz*dr;
    end
 end

% why are so few z values used?
% are the z and kz values without the bounds intended?