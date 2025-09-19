function [noise] = thermalnoise3D(Transf,Params,T,muchem)
    muchem = 0;
    
    %It is assumed T was given in nanokelvin
    Emax = 2*Params.kbol*T*1e-9/(Params.hbar*Params.w0);

    nxvals = 0:ceil(Emax*Params.w0/Params.wx-0.5);
    nyvals = 0:ceil(Emax*Params.w0/Params.wy-0.5);
    nzvals = 0:ceil(Emax*Params.w0/Params.wz-0.5);

    HOxEnergies = Params.wx*(nxvals+0.5)/Params.w0;
    HOyEnergies = Params.wy*(nyvals+0.5)/Params.w0;
    HOzEnergies = Params.wz*(nzvals+0.5)/Params.w0;

    [XEner,YEner,ZEner] = ndgrid(HOxEnergies,HOyEnergies,HOzEnergies);
    TotEner = XEner + YEner + ZEner;

    nBE = 1./(exp(2*(TotEner-muchem)/Emax)-1); % 2kBT = Emax, subtract chemical potential
    nBE(nBE<0) = 0; %negative occupanices are nonphysical

    r = normrnd(0,1,size(nBE));
    thet = rand(size(nBE));
    alpha = sqrt(nBE+1/2).*r.*exp(2*pi*1i*thet);
    
    %XX
%     NumTrials = 10000;
%     avg = 0;
%     for ii = 1:NumTrials
%         r = normrnd(0,1,size(nBE));
%         thet = rand(size(nBE));
%         alpha = sqrt(nBE+1/2).*r.*exp(2*pi*1i*thet);
%         avg = avg + sum(alpha(:).*conj(alpha(:)))/NumTrials;
%     end
%     disp(avg)
%     sum(nBE(:)+1/2)
    %XX

    alpha(TotEner>Emax)=0;

    [~,~,Z] = ndgrid(Transf.x,Transf.y,Transf.z);
    noise = zeros(size(Z));

    %Generate Hermites for 1D vector rather than 3D
    hermsx = []; hermsy = []; hermsz = [];
    xx = 1;
    for nx = nxvals
        hermsx(:,xx) = 1/sqrt(2^nx*factorial(nx)*sqrt(pi))*exp(-Transf.x.^2/2).*hermiteH(nx,Transf.x);
        xx = xx+1;
    end
    yy = 1;
    for ny = nyvals
        hermsy(:,yy) = 1/sqrt(2^ny*factorial(ny)*sqrt(pi))*exp(-Transf.y.^2/2).*hermiteH(ny,Transf.y);
        yy = yy+1;
    end
    zz = 1;
    for nz = nzvals
        hermsz(:,zz) = 1/sqrt(2^nz*factorial(nz)*sqrt(pi))*exp(-Transf.z.^2/2).*hermiteH(nz,Transf.z);
        zz = zz+1;
    end

    xx = 1;
    for nx = nxvals
        yy = 1;
        for ny = nyvals
            zz = 1;
            for nz = nzvals
                [HermX,HermY,HermZ] = ndgrid(hermsx(:,xx),hermsy(:,yy),hermsz(:,zz));
%                 Phi = alpha(xx,yy,zz)*HermX.*HermY.*hermsz(:,zz);
                Phi = alpha(xx,yy,zz)*HermX.*HermY.*HermZ;%*sqrt(Params.N);
                
                noise = noise + Phi;
                zz = zz + 1;
                % Norm = sum(abs(noise(:)).^2)*Transf.dx*Transf.dy*Transf.dz; % normalisation

            end
            yy = yy +1;            
        end
        xx = xx +1;
    end