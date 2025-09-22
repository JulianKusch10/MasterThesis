function [HTTable,densList] = calculate_Htherm(Params,Transf,VDk)
    addpath('lgwt'); %(File-Exch. #4540)
    addpath('GauLagWt'); %(File-Exch. #69136)

    Nk = 120; % Number of Gauss-Laguerre modes for k ∈ [k_min,k_max]
    Nu = 24; %Number of Gauss-Legendre modes for u ∈ [-1,1]

    Ndens = 1e3; %Number of density points to evaluate the integral at
    densList = logspace(-5, 6, Ndens);  %Params.N

    T = Params.T*1e-9/(Params.hbar*Params.w0); %assuming temperature given in nK
    beta   = 1/(Params.kbol*T);
    eps_dd = Params.add/Params.as;

    %Low momenta
    kxHO = 2*pi/(Params.ellx/Params.l0);
    kyHO = 2*pi/(Params.elly/Params.l0);
    kzHO = 2*pi/(Params.ellz/Params.l0);
    kLow = max([kxHO,kyHO,kzHO]);%0.017/(Params.add/Params.l0);
    kdd = 1/(Params.add/Params.l0);

    kIR = 0*0.5*kdd; % infrared cutoff k_min
    kUV = Inf; % ultraviolet cutoff k_max

    % Gauss-Legendre weights 
    [u,  wu]  = lgwt(Nu, -1, 1);                
   
    if isfinite(kUV) % Gauss-Legendre weights (finite UV cutoff) OR
        [k, wk] = lgwt(Nk, kIR, kUV);
    else  % Gauss-Laguerre weights (infinite UV cutoff)
        [k, wk]   = Gaulagwt(Nk);
        wk = wk.';
        k = k + kIR; % shift if kIR > 0
        wk = wk .* exp(k-kIR); % compensate the e^{-k} weight
    end
    k  = k(:);    wk = wk(:);
    u  = u(:).';  wu = wu(:).';

    %Kinetic energy
    tau_k  = k.^2/2;
            
    %We sum everything at once (every combination of Σ_k wk  Σ_u wu = Σ_kΣ_u wkwu
    [kMat,uMat] = ndgrid(k,u);
    tauMat = tau_k .* ones(1,Nu);
    
    P2 = 3*uMat.^2 - 1; % Legendre P_2(u)
    
    VkMat = Params.gs .* (1 + eps_dd*P2); %(depends on angle)

    % outer product for all weights
    % Wmat   = wk * wu.';
    Wmat   = wk .* wu;

    % -- allocate output 
    HTTable = zeros(1, Ndens);
    tic
    for nn = 1:Ndens %loop over densities
        rho = densList(nn);
    
        epsMat = sqrt(tauMat.*(tauMat + 2*rho*VkMat)); % Bogoliubov spectrum
        BE = 1 ./ (exp(beta*epsMat) - 1); % Bose factor
    
        integrand = kMat.^2.*VkMat .* tauMat ./ epsMat .* BE; % Nk × Nu
    
        Wmat(abs(imag(epsMat))>1e-15) = 0;

        HTTable(nn) = sum( integrand(:) .* Wmat(:) );  % double sum
    end
    toc
    
    HTTable = HTTable/(2*pi)^2; % put prefactor outside (one factor of (2pi) gone from phi integral)

    %%
    % figure(1)
    plot(log10(densList),log10(real(HTTable)),'LineStyle','-','Marker','o','Color','b','MarkerFaceColor','b','LineWidth',1,'MarkerSize',2,'DisplayName','Spherical');
    hold on;
    