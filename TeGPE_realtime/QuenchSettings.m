function Quench = QuenchSettings(Params,muchem_initial,muchem_final)

Quench.dt = 0.005; %Quench dt time in ms
Quench.dt = Quench.dt*Params.w0/1000; %Quench dt time in ms

Quench.T_initial = 10; %Temperature for the noise (in nanokelvin)
Quench.T_final = 10;

Quench.NumQuenches = 1;
Quench.QuenchIdxStart = 1;

Quench.eq_time = 0; %Equilibration time in ms
Quench.quench_time = 0; %Quench time in ms
Quench.hold_time = 1000; %Hold time in ms

Quench.as_initial = Params.as; %Initial scattering length in a0 
Quench.as_final = Params.as; %Final scattering length in a0

Quench.muchem_initial = muchem_initial;
Quench.muchem_final = muchem_final;

%save step
Quench.savestep = 400;

% --- Constructing quench vectors
eq_time = Quench.eq_time*Params.w0/1000;
quench_time = Quench.quench_time*Params.w0/1000;
hold_time = Quench.hold_time*Params.w0/1000;

tVec_Eq = 0:Quench.dt:eq_time;
tVec_Q = (eq_time+Quench.dt):Quench.dt:(eq_time+quench_time);
tVec_H = (eq_time+quench_time+Quench.dt):Quench.dt:(eq_time+quench_time+hold_time);

Quench.tSteps = length(tVec_Eq) + length(tVec_Q) + length(tVec_H);
Quench.tVec = [tVec_Eq,tVec_Q,tVec_H];

Quench.numsaves = floor(Quench.tSteps/Quench.savestep);

%temperature vector
Quench.Temp_vec = [Quench.T_initial*ones(1, length(tVec_Eq)), linspace(Quench.T_initial,Quench.T_final,length(tVec_Q)), Quench.T_final*ones(1, length(tVec_H))]; 

%linear as_quench
Quench.as_vec = [Quench.as_initial*ones(1, length(tVec_Eq)), linspace(Quench.as_initial,Quench.as_final,length(tVec_Q)), Quench.as_final*ones(1, length(tVec_H))]; 
Quench.gs_vec = 4*pi*Quench.as_vec/Params.l0;

%muchem vector
Quench.muchem_vec = [Quench.muchem_initial*ones(1, length(tVec_Eq)), linspace(Quench.muchem_initial,Quench.muchem_final,length(tVec_Q)), Quench.muchem_final*ones(1, length(tVec_H))]; 

% Quantum fluctuations need to be calculated over the quench
eps_dd_vec = Params.add./Quench.as_vec;
Q5Approx = 1;

if Q5Approx == 0
    yeps_vec = (1-eps_dd_vec)./(3*eps_dd_vec);
    Q5_vec = (3*eps_dd_vec).^(5/2).*( (8+26*yeps_vec+33*yeps_vec.^2).*sqrt(1+yeps_vec) + 15*yeps_vec.^3.*log((1+sqrt(1+yeps_vec))./sqrt(yeps_vec)) )/48;
    Q5_vec = real(Q5_vec);    
    Q5_vec(eps_dd_vec == 0) = 1;
    Q5_vec(eps_dd_vec == 1) = 3*sqrt(3)/2;
    Quench.gammaQF_vec = 128/3*sqrt(pi*(Quench.as_vec/Params.l0).^5).*Q5_vec;
else
    Q5_vec = 1+3*(eps_dd_vec).^2/2;
    Quench.gammaQF_vec = 128/3*sqrt(pi*(Quench.as_vec/Params.l0).^5).*Q5_vec;
end