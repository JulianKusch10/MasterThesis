function [psi,V,VDk] = Initialize(Params,Transf)

    format long
    X = Transf.X; Y = Transf.Y; Z = Transf.Z;    
    
    % == Potential == %
    V = 0.5*(Params.gx.*X.^2+Params.gy.*Y.^2+Params.gz*Z.^2); 

    % == Calulating the DDIs == %
    % [✓] Cylindrical (semianalytic) -- 1
	% [✓] Cylindrical infinite Z, polarized along x (analytic) -- 2
	% [✓] Cylindrical analytic, truncated -- 3
	% [✓] Spherical -- 4
	% Load -- 5
    % [✓] Square
	
	CutoffType = 6;

    switch CutoffType
		case 1 %Cylindrical (semianalytic)
            Zcutoff = Params.Lz/2;            
		    alph = acos((Transf.KX*sin(Params.theta)*cos(Params.phi)+Transf.KY*sin(Params.theta)*sin(Params.phi)+Transf.KZ*cos(Params.theta))./sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2));
		    alph(1) = pi/2;

            %Analytic part
            cossq = cos(alph).^2;
            VDk = cossq-1/3;
            sinsq = 1 - cossq;
            VDk = VDk + exp(-Zcutoff*sqrt(Transf.KX.^2+Transf.KY.^2)).*( sinsq .* cos(Zcutoff * Transf.KZ) - sqrt(sinsq.*cossq).*sin(Zcutoff * Transf.KZ) );
            
            %Nonanalytic part
            Params.Lr = 0.5*min(Params.Lx,Params.Ly);
		    Params.Nr = max(Params.Nx,Params.Ny);
		    [TransfRad] = setup_space_radial(Params); %morder really doesn't matter
		    VDkNon = VD_RadInt(TransfRad.kr,TransfRad.kz,TransfRad.Rmax,Zcutoff);
            
            %Interpolating the nonanalytic part onto 3D grid
            fullkr = [-flip(TransfRad.kr)',TransfRad.kr'];
		    [KR,KZ] = ndgrid(fullkr,TransfRad.kz);
		    [KX3D,KY3D,KZ3D] = ndgrid(ifftshift(Transf.kx),ifftshift(Transf.ky),ifftshift(Transf.kz));
		    KR3D = sqrt(KX3D.^2 + KY3D.^2);
		    fullVDK = [flip(VDkNon',2),VDkNon']';
		    VDkNon = interpn(KR,KZ,fullVDK,KR3D,KZ3D,'spline',0); %Last argument is -1/3 for full VDk. 0 for nonanalytic piece?
            VDkNon = fftshift(VDkNon);

            VDk = VDk + VDkNon;
            %
            VDk = 3*VDk; %finally fixing 1/3 inconsistency

		    save(sprintf('./Data/VDk_M.mat'),'VDk');
	        disp('Finished DDI')

		case 2 %Cylindrical infinite Z, polarized along x -- PRA 107, 033301 (2023)
		    alph = acos((Transf.KX*sin(Params.theta)*cos(Params.phi)+Transf.KY*sin(Params.theta)*sin(Params.phi)+Transf.KZ*cos(Params.theta))./sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2));
		    alph(1) = pi/2;
		    rhoc = max([abs(Transf.x),abs(Transf.y)]);
		    KR = sqrt(Transf.KX.^2+Transf.KY.^2);
		    func = @(n,u,v) v.^2./(u.^2+v.^2).*(v.*besselj(n,u).*besselk(n+1,v) - u.*besselj(n+1,u).*besselk(n,v));    
		    VDk = -0.5*func(0,KR*rhoc,abs(Transf.KZ)*rhoc) + (Transf.KX.^2./KR.^2 - 0.5).*func(2,KR*rhoc,abs(Transf.KZ)*rhoc);
		    VDk = (1/3)*(3*(cos(alph).^2)-1) - VDk;

            VDk(KR==0) = -1/3 + 1/2*abs(Transf.KZ(KR==0))*rhoc.*besselk(1,abs(Transf.KZ(KR==0))*rhoc);
		    VDk(Transf.KZ==0) = 1/6 + (Transf.KX(Transf.KZ==0).^2-Transf.KY(Transf.KZ==0).^2)./(KR(Transf.KZ==0).^2).*(1/2 - besselj(1,KR(Transf.KZ==0)*rhoc)./(KR(Transf.KZ==0)*rhoc));
		    VDk(1,1,1) = 1/6;
            
            VDk = 3*VDk; %finally fixing 1/3 inconsistency

        case 3 %Cylindrical Z, polarized along x [PRA 107, 033301 (2023)] -- then truncate along z direction
		    alph = acos((Transf.KX*sin(Params.theta)*cos(Params.phi)+Transf.KY*sin(Params.theta)*sin(Params.phi)+Transf.KZ*cos(Params.theta))./sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2));
		    alph(1) = pi/2;
		    rhoc = max([abs(Transf.x),abs(Transf.y)]);
		    KR = sqrt(Transf.KX.^2+Transf.KY.^2);
		    func = @(n,u,v) v.^2./(u.^2+v.^2).*(v.*besselj(n,u).*besselk(n+1,v) - u.*besselj(n+1,u).*besselk(n,v));    
		    VDk = -0.5*func(0,KR*rhoc,abs(Transf.KZ)*rhoc) + (Transf.KX.^2./KR.^2 - 0.5).*func(2,KR*rhoc,abs(Transf.KZ)*rhoc);
		    VDk = (1/3)*(3*(cos(alph).^2)-1) - VDk;

            VDk(KR==0) = -1/3 + 1/2*abs(Transf.KZ(KR==0))*rhoc.*besselk(1,abs(Transf.KZ(KR==0))*rhoc);
		    VDk(Transf.KZ==0) = 1/6 + (Transf.KX(Transf.KZ==0).^2-Transf.KY(Transf.KZ==0).^2)./(KR(Transf.KZ==0).^2).*(1/2 - besselj(1,KR(Transf.KZ==0)*rhoc)./(KR(Transf.KZ==0)*rhoc));
		    VDk(1,1,1) = 1/6;
            
            %z cutoff
            Zcutoff = Params.Lz/2;
            VDr = fftshift(fftn(VDk));
            VDr(abs(Transf.Z)>Zcutoff) = 0;
            VDk = ifftn(ifftshift(VDr));
            
            VDk = 3*VDk; %finally fixing 1/3 inconsistency

        case 4 %Spherical
            Rcut = min(Params.Lx/2,Params.Ly/2,Params.Lz/2);
            alph = acos((Transf.KX*sin(Params.theta)*cos(Params.phi)+Transf.KY*sin(Params.theta)*sin(Params.phi)+Transf.KZ*cos(Params.theta))./sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2));
		    alph(1) = pi/2;

            K = sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2);
            VDk = (cos(alph).^2-1/3).*(1 + 3*cos(Rcut*K)./(Rcut^2.*K.^2) - 3*sin(Rcut*K)./(Rcut^3.*K.^3));

            VDk = 3*VDk; %finally fixing 1/3 inconsistency

        case 5
            VDk = load(sprintf('./Data/VDk_M.mat'));
		    VDk = VDk.VDk;

        case 6 %Rectangular
            alph = acos((Transf.KX*sin(Params.theta)*cos(Params.phi)+Transf.KY*sin(Params.theta)*sin(Params.phi)+Transf.KZ*cos(Params.theta))./sqrt(Transf.KX.^2+Transf.KY.^2+Transf.KZ.^2));
		    alph(1) = pi/2;
            VDk = (cos(alph).^2-1/3);
            VDk(1,1,1) = 0;

            %xyz cutoff
            Xcutoff = Params.Lx/2;
            Ycutoff = Params.Ly/2;
            Zcutoff = Params.Lz/2;
            VDr = fftshift(fftn(VDk));
            VDr(abs(Transf.X)>Xcutoff) = 0;
            VDr(abs(Transf.Y)>Ycutoff) = 0;
            VDr(abs(Transf.Z)>Zcutoff) = 0;
            VDk = ifftn(ifftshift(VDr));
            VDk = 3*VDk; %finally fixing 1/3 inconsistency
            
        otherwise
		    disp('Choose a valid DDI!')
            return
    end


    % == Setting up the initial wavefunction == %
    loadstate = 0;
    
    if loadstate == 1
        loadnumber = 1250;
    %     load(sprintf('./Data/Seed/psi_%i.mat',loadnumber),'psi');
        load(sprintf('./Data/Run_004/psi_gs.mat'),'psi');
    
        Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
        psi = sqrt(Params.N)*psi/sqrt(Norm);
    
    else
        ellx = sqrt(Params.hbar/(Params.m*Params.wx))/Params.l0; 
        elly = sqrt(Params.hbar/(Params.m*Params.wy))/Params.l0;
        ellz = sqrt(Params.hbar/(Params.m*Params.wz))/Params.l0; 
    
        Rx = 8*ellx; Ry = 8*elly; Rz = 8*ellz;
        X0 = 0.0*Transf.Xmax; Y0 = 0.0*Transf.Ymax; Z0 = 0*Transf.Zmax;
    
        psi = exp(-(X-X0).^2/Rx^2-(Y-Y0).^2/Ry^2-(Z-Z0).^2/Rz^2);
        cur_norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
        psi = psi/sqrt(cur_norm);
        
        % psi = zeros(size(psi)); %xxxxxx

        %add some noise
        r = normrnd(0,1,size(X));
        theta = rand(size(X));
        noise = r;%.*exp(2*pi*1i*theta);
        psi = psi + 0.005*noise;
    
        Norm = sum(abs(psi(:)).^2)*Transf.dx*Transf.dy*Transf.dz;
        psi = sqrt(Params.N)*psi/sqrt(Norm);
    end