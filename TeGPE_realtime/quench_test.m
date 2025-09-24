function quench_test()

    Quench.dt = 0.005;
    Quench.eq_time = 100; %Equilibration time i
    Quench.quench_time = 100; %Quench time in m
    Quench.hold_time = 1000; %Hold time in ms

    eq_time = Quench.eq_time;
    quench_time = Quench.quench_time;
    hold_time = Quench.hold_time;


    tVec_Eq = 0:Quench.dt:eq_time;
    tVec_Q = (eq_time+Quench.dt):Quench.dt:(eq_time+quench_time);
    tVec_H = (eq_time+quench_time+Quench.dt):Quench.dt:(eq_time+quench_time+hold_time);


    Quench.tSteps = length(tVec_Eq) + length(tVec_Q) + length(tVec_H);
    Quench.tVec = [tVec_Eq,tVec_Q,tVec_H];

    Quench.as_initial = 130;
    Quench.as_final = 130; 
    Quench.as_vec = [Quench.as_initial*ones(1, length(tVec_Eq)), linspace(Quench.as_initial,Quench.as_final,length(tVec_Q)), Quench.as_final*ones(1, length(tVec_H))];
    plot(Quench.tVec, Quench.as_vec);
    ylim([0, 150])
end 