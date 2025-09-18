function Quench = QuenchSettings(Params, taskid)
    Quench.NumTrials = 30; %Number of independent trials

    Quench.quench_time = 250; %Quench time in ms

    Quench.as = Params.as; %Quench scattering length in a0
    
    Quench.addnoise = 1;
    Quench.Temp = 50; %Temperature in nanoKelvin
    
    % --- Save steps
    Quench.savestep = 100;

    % --- Constructing quench vectors
    quench_time = Quench.quench_time*Params.w0/1000;
    tVec_Q = 0:Params.dt:quench_time;
    
    Quench.tSteps = length(tVec_Q);
    Quench.as_vec = Quench.as*ones(1, length(tVec_Q));
    Quench.gs_vec = 4*pi*Quench.as_vec/Params.l0;
    Quench.tVec = tVec_Q;
    
    Quench.ObservSteps = floor(Quench.tSteps/Quench.savestep) + 2;
    % Quantum fluctuations need to be calculated over the quench
    eps_dd_vec = Params.add./Quench.as_vec;
    yeps_vec = (1-eps_dd_vec)./(3*eps_dd_vec);
    Q5_vec = (3*eps_dd_vec).^(5/2).*( (8+26*yeps_vec+33*yeps_vec.^2).*sqrt(1+yeps_vec) + 15*yeps_vec.^3.*log((1+sqrt(1+yeps_vec))./sqrt(yeps_vec)) )/48;
    Q5_vec = real(Q5_vec);
    Q5_vec(eps_dd_vec == 0) = 1;
    Q5_vec(eps_dd_vec == 1) = 3*sqrt(3)/2;
    Quench.gammaQF_vec = 128/3*sqrt(pi*(Quench.as_vec/Params.l0).^5).*Q5_vec;    
%
