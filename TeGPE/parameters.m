function [Params] = parameters(aidx,asList)
%%--%% Parameters %%--%%
%========= Simulation =========%
pert=0; %0 = no perturbation during real-time, 1=perturbation
%method=1; %0 = normal dipolar potential, 1=spherical cut-off, 2=cylindrical cut-off

%Energy tolerance
Params.Etol=5e-10;
Params.rtol = 1e-6;
Params.theta = 0; %theta=0 dipoles along z ,pi/2 dipoles along x,
Params.phi= 0; %azimuthal angle
Params.cut_off=3e7; %sometimes the imaginary time gets a little stuck
              %even though the solution is good, this just stops it going on forever

%========= Constants =========%
Params.hbar = 1.0545718e-34; %Planck constant [J.s]
Params.kbol = 1.38064852e-23; %Boltzmann Constant [J/K]
Params.mu0 = 1.25663706212e-6; %Vacuum Permeability [N/A^2]  --
Params.muB = 9.274009994e-24; %Bohr Magneton [J/T]
Params.a0 = 5.2917721067e-11; %Bohr radius [m]
Params.m0 = 1.660539066e-27; %Atomic mass [kg]
Params.w0 = 2*pi*61.63158647; %angular frequency unit [s^-1]
%62.39247025-> Dysprosium 162
%61.63158647-> Dysprosium 164
%60.88903723-> Erbium 166
Params.mu0factor = 0.3049584233607396;% =(m0/me)*pi*alpha^2 -- me=mass of electron, alpha=fine struct. const.
                               % mu0=mu0factor *hbar^2*a0/(m0*muB^2)
%=============================%

%Project 
Params.name = "test2";
Params.save_plots = 1;
Params.save_params = 1;

%Stopping 
Params.stop_relres_flag = 1;
Params.stop_relres = 10^(-3);


%Number of points in each direction
Params.Nx = 16;
Params.Ny = 16;
Params.Nz = 16;

%Dimensions (in units of l0)
Params.Lx = 24;
Params.Ly = 46;
Params.Lz = 14;

%Masses
Params.m = 164*Params.m0;
Params.l0 = sqrt(Params.hbar/(Params.m*Params.w0)); %Defining a harmonic oscillator length

%Atom numbers
% Params.ppum = 2500; %particles per micron
% Params.N = Params.Lz*Params.ppum*l0*1e6;

Params.N = 8*10^4;%2.4*10^4;

%Dipole lengths (units of muB)
Params.mu = 9.93*Params.muB;

%scattering lengths
 %asList(taskid)
Params.as = 90*Params.a0;

%trapping frequencies
Params.wx = 2*pi*50;
Params.wy = 2*pi*20;
Params.wz = 2*pi*150;

%Time step
Params.dt = 1e-3;
Params.mindt = 2e-8; %Minimum size for a time step using adaptive dt

%Thermal GPE
Params.T = 95; %Temperature in nanokelvin

%================ Parameters defined by those above ================%
%Contact interaction strength (units of l0/m)
Params.gs = 4*pi*Params.as/Params.l0;

%Dipole lengths
Params.add = Params.mu0*Params.mu^2*Params.m/(12*pi*Params.hbar^2);

% == Calculating quantum fluctuations == %
eps_dd = Params.add/Params.as;
Q5Approx = 1;

if Q5Approx == 1
    Q5 = 1+(3/2)*eps_dd^2;
else
    if eps_dd == 0
        Q5 = 1;
    elseif eps_dd == 1
        Q5 = 3*sqrt(3)/2;
    else
        yeps = (1-eps_dd)/(3*eps_dd);
        Q5 = (3*eps_dd)^(5/2)*( (8+26*yeps+33*yeps^2)*sqrt(1+yeps) + 15*yeps^3*log((1+sqrt(1+yeps))/sqrt(yeps)) )/48;
        Q5 = real(Q5);
    end
end

Params.gammaQF = 128/3*sqrt(pi*(Params.as/Params.l0)^5)*Q5;
if Params.as == 0
    Params.gammaQF = 0;
end

%DDI strength
Params.gdd = 4*pi*Params.add/Params.l0; %sometimes the 12 is a 4 --> depends on how Vdk (DDI) is defined
                                  % usually Vdk = gdd(3*cos^2theta-1)
                                  % but we often write Vdk=(cos^2theta-1/3)
                                  % so a factor of 3 is absorbed into gdd
%Trap lengths
Params.ellx = sqrt(Params.hbar/(Params.m*Params.wx));
Params.elly = sqrt(Params.hbar/(Params.m*Params.wy));
Params.ellz = sqrt(Params.hbar/(Params.m*Params.wz));

%Trap gamma
Params.gx=(Params.wx/Params.w0)^2;
Params.gy=(Params.wy/Params.w0)^2;
Params.gz=(Params.wz/Params.w0)^2;

%Loading the rest into Params
% Params.hbar = hbar; Params.kbol = kbol; Params.mu0 = mu0; Params.muB = muB; Params.a0 = a0;
% Params.w0 = w0; Params.l0 = l0;