function [Params] = parameters(taskid)
%%--%% Parameters %%--%%
%========= Simulation =========%
pert=0; %0 = no perturbation during real-time, 1=perturbation
%method=1; %0 = normal dipolar potential, 1=spherical cut-off, 2=cylindrical cut-off

%Energy tolerance
Params.Etol=5e-10;
Params.rtol = 1e-8;
Params.cut_off=3e7; %sometimes the imaginary time gets a little stuck
              %even though the solution is good, this just stops it going on forever

%========= Constants =========%
hbar = 1.0545718e-34; %Planck constant [J.s]
kbol = 1.38064852e-23; %Boltzmann Constant [J/K]
mu0 = 1.25663706212e-6; %Vacuum Permeability [N/A^2]  --
muB = 9.274009994e-24; %Bohr Magneton [J/T]
a0 = 5.2917721067e-11; %Bohr radius [m]
m0 = 1.660539066e-27; %Atomic mass [kg]
w0 = 2*pi*61.63158647; %angular frequency unit [s^-1]
%62.39247025-> Dysprosium 162
%61.63158647-> Dysprosium 164
%60.88903723-> Erbium 166
mu0factor = 0.3049584233607396;% =(m0/me)*pi*alpha^2 -- me=mass of electron, alpha=fine struct. const.
                               % mu0=mu0factor *hbar^2*a0/(m0*muB^2)
%=============================%

%Number of points in each direction
Params.Nx = 16;
Params.Ny = 16;
Params.Nz = 16;

%Dimensions (in units of l0)
Params.Lx = 30; %24;
Params.Ly = 20; %40;
Params.Lz = 15; %14;

%Masses
Params.m = 164*m0;
l0 = sqrt(hbar/(Params.m*w0)); %Defining a harmonic oscillator length

Params.N = 6.3*10^4; %10*10^4;


%Dipole lengths (units of muB)
Params.mu = 9.93*muB;

%scattering lengths
Params.as = 88*a0; %86*a0; %aslist(taskid)

%trapping frequencies
Params.wx = 2*pi*33; %2*pi*47;
Params.wy = 2*pi*84.6;%2*pi*21;
Params.wz = 2*pi*167;%2*pi*150;

% dipole angles
angles_deg = 0:5:50;
Params.theta = 0;%pi*angles_deg(taskid)/180; %theta=0 dipoles along z, pi/2 dipoles along x,
Params.phi = 0; %azimuthal angle


%Time step
Params.dt = 1e-3;
Params.mindt = 1e-6; %Minimum size for a time step using adaptive dt

%================ Parameters defined by those above ================%
%Contact interaction strength (units of l0/m)
Params.gs = 4*pi*Params.as/l0;

%Dipole lengths
Params.add = mu0*Params.mu^2*Params.m/(12*pi*hbar^2);

% == Calculating quantum fluctuations == %
eps_dd = Params.add/Params.as; 
if eps_dd == 0
    Q5 = 1;
elseif eps_dd == 1
    Q5 = 3*sqrt(3)/2;
else
    yeps = (1-eps_dd)/(3*eps_dd);
    Q5 = (3*eps_dd)^(5/2)*( (8+26*yeps+33*yeps^2)*sqrt(1+yeps) + 15*yeps^3*log((1+sqrt(1+yeps))/sqrt(yeps)) )/48;
    Q5 = real(Q5);
end

Params.gammaQF = 128/3*sqrt(pi*(Params.as/l0)^5)*Q5;
if Params.as == 0
    Params.gammaQF = 0;
end

%DDI strength
Params.gdd = 12*pi*Params.add/l0; %sometimes the 12 is a 4? --> depends on how Vdk (DDI) is defined

%Trap gamma
Params.gx=(Params.wx/w0)^2;
Params.gy=(Params.wy/w0)^2;
Params.gz=(Params.wz/w0)^2;

%Loading the rest into Params
Params.hbar = hbar; Params.kbol = kbol; Params.mu0 = mu0; Params.muB = muB; Params.a0 = a0;
Params.w0 = w0; Params.l0 = l0;