function HThermFit = calculate_Htherm(Params,kCutoff)

    if Params.T == 0 %skips the calculation if temperature is zero, but still works with the rest of code
        HThermFit.HTFitModel = @(~, xdata) zeros(size(xdata));
    
        zfun = @(~,~) 0;       % returns 0 for any (T, a_s)
        coeffFits.A = zfun;    coeffFits.B = zfun;    coeffFits.C = zfun;
    
        HThermFit.coeffFit      = coeffFits;
        HThermFit.coeffFitType  = [];
        return
    end

    %%
    asList = [132:4:180]*Params.a0; %Scattering length in a0
    TList = [5:4:120]; %Temperatures in nanoKelvin
    gsList = 4*pi*asList/Params.l0;

    [TT,AS] = ndgrid(TList,asList);

    lsqoptions = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000);
     
    InitGuess = [2.0 0.5];
    LB = [0 0];
    UB = [8.0 4.0];

    %Guesses for other models (see below)
    % InitGuess = [2.0 0.5 0.5]; LB = [0 0 0]; UB = [9.0 4.0 2.0];
    % InitGuess = [2.0 0.5 0.5 0.2]; LB = [0 0 0 0]; UB = [9.0 4.0 2.0 1.0];
    % InitGuess = [2.0 0.5 0.5 0.2 0.05]; LB = [0 0 0 0 0]; UB = [9.0 4.0 2.0 1.0 1.0];

    NumParams = length(InitGuess);
    ParamList = zeros(NumParams,length(TList),length(asList));

    for ii = 1:length(TList)
        
        Params.T = TList(ii);
        for jj = 1:length(asList)
            Params.as = asList(jj);
            Params.gs = gsList(jj);
    
            [HTTable,densList] = cylindrical_Htherm(Params,180,kCutoff);

            HTFitModel = @(x,xdata) x(1)*exp(-x(2)*xdata.^0.5);
            %Other models:
            % HTFitModel = @(x,xdata) x(1)*exp(-x(2)*xdata.^0.5-x(3)*xdata);
            % HTFitModel = @(x,xdata) x(1)./(1 + x(2)*xdata.^0.5);
            % HTFitModel = @(x,xdata) x(1)./(1 + x(2)*xdata.^0.5 + x(3)*xdata + x(4)*xdata.^1.5);
            % HTFitModel = @(x,xdata) x(1)./(1 + x(2)*xdata.^0.5 + x(3)*xdata + x(4)*xdata.^1.5 + x(5)*xdata.^2.0);

            if ii > 1
                if jj > 1
                    InitGuess = ParamList(:,ii,jj-1);
                else
                    InitGuess = ParamList(:,ii-1,jj);
                end
            end

            xfit = lsqcurvefit(HTFitModel, InitGuess, densList', HTTable',LB, UB, lsqoptions);
            ParamList(:,ii,jj) = xfit;
        end
    end
    
    asList_extended = linspace(90,220,100);
    TList_extended = linspace(2,150,100);

    [TT_E,AS_E] = ndgrid(TList_extended,asList_extended);

    coeffFits = struct(); %initialize an empty structure
    fitLabels = {'A','B','C','D','E'};

    for ii = 1:NumParams
        Z = squeeze(ParamList(ii,:,:));

        % Flatten for fitting
        xData = AS(:)/Params.a0;
        yData = TT(:);
        zData = Z(:);

        if ii == 1
            XY = [xData, yData];
            zVec = zData;

            % poly22
            modelFun = @(p, XY) p(1) + p(2)*XY(:,1) + p(3)*XY(:,2) + p(4)*XY(:,1).^2 + p(5)*XY(:,1).*XY(:,2) + p(6)*XY(:,2).^2;
            initGuess = zeros(1,6);
            opts = optimoptions('lsqcurvefit', 'Display', 'off');
            pFit = lsqcurvefit(modelFun, initGuess, XY, zVec, [], [], opts);

            % Wrap into function handle
            coeffFit = @(x, y) reshape(modelFun(pFit, [x(:), y(:)]), size(x));
            coeffFunc = modelFun;

            % Store parameter info
            coeffFitType = struct('model', 'poly22', 'params', pFit);
        elseif ii>=2
            % Prepare input data for lsqcurvefit
            XY = [xData, yData];       % N x 2 matrix
            zVec = zData;              % Flattened z values (already done)
        
            % Define model: a / (b*y + c*sqrt(y) + d) + e*exp(-f*x)
            modelFun = @(p, XY) p(1) ./ (p(2)*XY(:,2) + p(3)*sqrt(XY(:,2)) + p(4)) + p(5)*exp(-p(6)*XY(:,1));  % p = [a, b, c, d, e, f]
        
            % Starting point and bounds
            startVals = [0.1, 0.05, 0.05, 0, 0.05, 0.01];
            lb = [0, 0, 0, 0, 0, 0];  % All parameters ≥ 0
            ub = [];
        
            opts = optimoptions('lsqcurvefit', 'Display', 'off');
            pFit = lsqcurvefit(modelFun, startVals, XY, zVec, lb, ub, opts);
        
            % Wrap fit into anonymous function for consistency with `fit` output
            coeffFit = @(x, y) reshape(modelFun(pFit, [x(:), y(:)]),size(x));
            coeffFunc = modelFun;
            % Store parameter vector in structure for tracking
            coeffFitType = struct('model', 'custom_rational_exp', 'params', pFit);
        end
        
        SSR = abs(coeffFit(AS/Params.a0, TT)-Z).^2; %Maybe improve this with R^2 instead: R^2=1-
        % SSR = sum(SSR(:))

        % Plot
        Zfit_extended = coeffFit(AS_E, TT_E);
        subplot(2,3,ii)
        surf(AS/Params.a0, TT, Z,'FaceColor','r','EdgeAlpha',0.6,'DisplayName',fitLabels{ii})
        hold on
        surf(AS_E, TT_E, Zfit_extended,'FaceAlpha',0.6,'EdgeAlpha',0.6)

        xlabel('Scattering length -- x')
        ylabel('Temperature -- y')
        
        % Store in structure with meaningful field names
        coeffFits.(fitLabels{ii}) = coeffFit;
        coeffFitTypes.(fitLabels{ii}) = coeffFitType;
        coeffFitFunc.(fitLabels{ii}) = coeffFunc;
    end
    % legend
    HThermFit.HTFitModel = HTFitModel;
    HThermFit.coeffFit = coeffFits;
    HThermFit.coeffFitType = coeffFitTypes;
    HThermFit.coeffFitFunc = coeffFunc;
end

%% Function that actually calculates HTherm using
function [HTTable,densList] = cylindrical_Htherm(Params,NLaguerre,kCutoff)
    addpath('lgwt'); %(File-Exch. #4540)
    addpath('GauLagWt'); %(File-Exch. #69136)

    NLegendre = 24; %Number of Gauss-Legendre modes for u ∈ [-1,1]

    Ndens = 2e3; %Number of density points to evaluate the integral at
    densList = logspace(-5.01, 4.01, Ndens); %logspace(-5, 6, Ndens);  %Params.N

    T = Params.T*1e-9/(Params.hbar*Params.w0); %assuming temperature given in nK
    beta   = 1/(Params.kbol*T);
    eps_dd = Params.add/Params.as;

    %Low momenta
    kxHO = 2*pi/(Params.ellx/Params.l0);
    kyHO = 2*pi/(Params.elly/Params.l0);
    kzHO = 2*pi/(Params.ellz/Params.l0);
    kdd = 1/(Params.add/Params.l0);

    %cutoffs
    switch kCutoff
        case 0
            kRhoIR = 0*sqrt(kxHO^2+kyHO^2);
            kRhoUV = Inf;
            kzIR = 0*kzHO;
            kzUV = Inf;
        case 1
            kRhoIR = 0.7*sqrt(kxHO^2+kyHO^2);
            kRhoUV = Inf;
            kzIR = 0.7*kzHO;
            kzUV = Inf;
    end

    %Generate kRho weights
    if isfinite(kRhoUV) % Gauss-Legendre weights (finite UV cutoff) OR
        [kRho, wkRho] = lgwt(NLegendre, kRhoIR, kRhoUV);
    else  % Gauss-Laguerre weights (infinite UV cutoff)
        [kRho, wkRho] = Gaulagwt(NLaguerre);
        wkRho = wkRho.';
        kRho = kRho + kRhoIR; % shift if kIR > 0
        wkRho = wkRho .* exp(kRho-kRhoIR); % compensate the e^{-k} weight
    end

    %Generate kz weights
    if isfinite(kzUV) % Gauss-Legendre weights (finite UV cutoff) OR
        [kz, wkz] = lgwt(NLegendre, kzIR, kzUV);
    else  % Gauss-Laguerre weights (infinite UV cutoff)
        [kz, wkz]   = Gaulagwt(NLaguerre);
        wkz = wkz.';
        kz = kz + kzIR; % shift if kIR > 0
        wkz = wkz .* exp(kz-kzIR); % compensate the e^{-k} weight
    end

    % Including Jacobian (kRho) into kRho weights
    wkRho = 2 * wkRho .* kRho; %factor of 2 because kz really goes from (-Inf,Inf) but we integrate [0,Inf)

    % Constructing terms
    [kRhoMat, kzMat] = ndgrid(kRho, kz);
    tauMat = (kRhoMat.^2 + kzMat.^2) / 2; %Kinetic energy
    cos2   = kzMat.^2 ./ (kRhoMat.^2 + kzMat.^2);
    P2Mat  = 3*cos2 - 1; %Legendre P2
    VkMat  = Params.gs .* (1 + eps_dd*P2Mat);  % interaction Vk

    %outer product for all weights
    Wmat   = wkRho * wkz.';
    
    % HTherm table to be interpolated
    HTTable = zeros(1, Ndens);
        
    for nn = 1:Ndens
        n   = densList(nn);
        epsMat = sqrt( tauMat .* ( tauMat + 2*n*VkMat ) );
        BE     = 1 ./ ( exp(beta*epsMat) - 1 );
    
        integrand = VkMat .* tauMat ./ epsMat .* BE;   % Nk_p × Nk_z
        
        Wmat(abs(imag(epsMat))>1e-12) = 0;
        % integrand = real(integrand);

        HTTable(nn) = sum( integrand(:) .* Wmat(:) );  % double sum

        % igcomplex = zeros(size(integrand));
        % igcomplex(abs(imag(epsMat))>1e-10) = 1;
        % subplot(4,4,nn)
        % pcolor(kzMat/kdd,kRhoMat/kdd,igcomplex','EdgeColor','none')
        % xlabel('$k_\rho$')
        % ylabel('$k_z$')
        % colormap('nebula')
        % text(0.025,0.015,sprintf('$n=%.2f \\mu\\mathrm{m}^{-3}$',n))
        % xlim([0,0.1])
        % ylim([0,0.02])
        % clim([0,1])
    end
    HTTable = HTTable / (2*pi)^2;    % one factor of (2pi) gone from phi integral

end

    