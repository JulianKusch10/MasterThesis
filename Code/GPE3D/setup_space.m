function [Transf] = setup_space(Params)
Transf.Xmax = 0.5*Params.Lx;
Transf.Ymax = 0.5*Params.Ly;
Transf.Zmax = 0.5*Params.Lz;

Nz = Params.Nz; Nx = Params.Nx; Ny = Params.Ny;

% Fourier grids
x = linspace(-0.5*Params.Lx,0.5*Params.Lx,Params.Nx); %Symmetric grid!!
Kmax = pi*Params.Nx/Params.Lx;
kx = linspace(-Kmax,Kmax,Nx+1);
kx = kx(1:end-1); dkx = kx(2)-kx(1);
kx = fftshift(kx);

y = linspace(-0.5*Params.Ly,0.5*Params.Ly,Params.Ny);
Kmax = pi*Params.Ny/Params.Ly;
ky = linspace(-Kmax,Kmax,Ny+1);
ky = ky(1:end-1); dky = ky(2)-ky(1);
ky = fftshift(ky);

z = linspace(-0.5*Params.Lz,0.5*Params.Lz,Params.Nz);
Kmax = pi*Params.Nz/Params.Lz;
kz = linspace(-Kmax,Kmax,Nz+1);
kz = kz(1:end-1); dkz = kz(2)-kz(1);
kz = fftshift(kz);

[Transf.X,Transf.Y,Transf.Z]=ndgrid(x,y,z);
[Transf.KX,Transf.KY,Transf.KZ]=ndgrid(kx,ky,kz);
Transf.x = x; Transf.y = y; Transf.z = z;
Transf.kx = kx; Transf.ky = ky; Transf.kz = kz;
Transf.dx = x(2)-x(1); Transf.dy = y(2)-y(1); Transf.dz = z(2)-z(1);
Transf.dkx = dkx; Transf.dky = dky; Transf.dkz = dkz;