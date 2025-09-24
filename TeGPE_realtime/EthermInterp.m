function Ftherm = EthermInterp(Params,HThermFit)        

    if Params.T == 0
        Ftherm = @(n) zeros(size(n));
        return
    end


    %% Calculate the antiderivative
    % In order to calculate the thermal energy, we also require the
    % antiderivative of our function w.r.t. n

    T = Params.T;
    as = Params.as/Params.a0;
        
    numcoeff = length(fieldnames(HThermFit.coeffFit));
    clabels = {'A','B','C','D','E'};
    
    coeffs = [];
    for ii = 1:numcoeff
        cmodel = HThermFit.coeffFit.(clabels{ii});
        coeffs(ii) = cmodel(as,T);
    end

    logdenslist = linspace(-9,6,1e3); %We can use a very broad density list now
    densList = 10.^logdenslist; 
    logdenslist = [coeffs(1), logdenslist]; %Because A is the zero density limit
    densList = [0, densList];
    HT = HThermFit.HTFitModel(coeffs,densList);

    EthTable = cumtrapz(densList,HT);
    Ftherm = griddedInterpolant(densList, EthTable, 'pchip', 'linear'); %F_T