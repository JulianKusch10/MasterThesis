function HT = HTherm(HThermFit,Params,psi)
    
    T = Params.T;
    as = Params.as/Params.a0;
    
    n = abs(psi).^2;
    
    numcoeff = length(fieldnames(HThermFit.coeffFit));
    clabels = {'A','B','C','D','E'};
    
    coeffs = [];
    for ii = 1:numcoeff
        cmodel = HThermFit.coeffFit.(clabels{ii});
        coeffs(ii) = cmodel(as,T);
    end
    
    HT = HThermFit.HTFitModel(coeffs,n);
    HT(HT<0) = 0; %Eliminate unphysical negative thermal energy from the fit function