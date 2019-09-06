clear
N = 5;
K = 32;

useSBP = 1;

alpha = 0; % 0 = upwinding

% plotting points
Npts = 100; a = -K/2; b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';

% mass quadrature
rq = []; wq = [];
h = 1; 
for e = 1:K    
    % full quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];    
end

h = 2/K;

VN = GDVDM(N,K,rq);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

% possibly reduced quadrature versions
Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);

VX = VX*h;
M = h*M;
% Q = Q;

% plot(VX,(M\Q)*exp(VX),'o')
% hold on
% plot(VX,exp(VX),'--')

M2 = M;
pskip = N+2;
inds = pskip:size(M,2)-(pskip)+1;
MDiag = diag(sum(M2,2));
M2(inds,:) = MDiag(inds,:);

sM2 = sum(M2,2);
M2(inds,inds) = diag(sM2(inds));
% M2 = diag(sum(M,2));

if useSBP    
    [MSBP, QSBP, X] = CSBPp3(K+1);
end
% M = MSBP; 
% Q = QSBP;
% M2 = MSBP;

%% interior dissipation

if 0
    mid = round(size(Q,2)/2);
    aa = Q(mid,mid+(1:N)); %/ M2(mid,mid);
    
    aaSBP = QSBP(mid,mid+(1:N));% / MSBP(mid,mid);
    Npts = 100;
    Lmin = 1/Npts;
    Lmax = pi-1/Npts;
    Lvec = linspace(Lmin,Lmax,Npts);
    for i = 1:Npts
        ell = Lvec(i);
        ilamhist_interior(i) = 2*sum(aa.*sin(ell*(1:N)));
        ilamhist_interior_SBP(i) = 2*sum(aaSBP.*sin(ell*(1:N)));
    end
    
    figure(1)
    loglog(Lvec,abs(ilamhist_interior-Lvec),'bo-','linewidth',2,'DisplayName','Lumped GD');
    hold on
    loglog(Lvec,abs(ilamhist_interior_SBP-Lvec),'rs-','linewidth',2,'DisplayName','SBP');
    % plot(Lvec,Lvec,'k--')
    loglog(Lvec,Lvec.^(N+2),'k--','DisplayName','h^{p+2}')
    loglog(Lvec,Lvec.^(2*N+1),'k--','DisplayName','h^{2p+1}')
    
    title('Interior node dispersion')
    legend show
    return
end

%%

e0 = zeros(K+1,1); e0(1) = 1;
eN = zeros(K+1,1); eN(K+1) = 1;

Npts = 50;
Lmin = 1/Npts;
Lmax = pi-1/Npts;
Lvec = linspace(Lmin,Lmax,Npts);
for i = 1:length(Lvec)
    
    ell = Lvec(i);
    k = 1i*ell*(K+1); % length of domain
    S = -alpha/2*eN*(eN' - exp(k)*e0') + (2-alpha)/2*e0*(e0' - exp(-k)*eN');
    
    % ---------    
    [~, D] = eig(Q+S,1i*M/2);
    lam = diag(D) / size(S,2);
    rlam = real(lam);
    ilam = imag(lam);
    
    if i < 3
        [~,id] = min(abs(rlam-ell));
    else % find closest linear extrapolation assuming dL constant
        newlam = rlamhist(i-2) + (rlamhist(i-1) - rlamhist(i-2))*2;
        inewlam = ilamhist(i-2) + (ilamhist(i-1) - ilamhist(i-2))*2;
        [~,id] = min(abs(rlam-newlam).^2 + abs(ilam-inewlam).^2); 
    end
    rlamhist(i) = rlam(id);    
    ilamhist(i) = ilam(id);    
    
    % -------------
    
    I = eye(size(M));
    Ar = M2\(M2-M);
%     invMlump = ((I+Ar+Ar.^2+Ar.^3)/M2);
%     Mlump = inv(invMlump);
    Mlump = M2;

    [~, D] = eig(Q+S,1i*Mlump/2);
    lam = diag(D) / size(S,2);
    rlam = real(lam);
    ilam = imag(lam);
    
    if i < 3
        [~,id] = min(abs(rlam-ell));
    else % find closest linear extrapolation assuming dL constant
        newlam = rlamhist2(i-2) + (rlamhist2(i-1) - rlamhist2(i-2))*2;
        inewlam = ilamhist2(i-2) + (ilamhist2(i-1) - ilamhist2(i-2))*2;
        [~,id] = min(abs(rlam-newlam).^2 + abs(ilam-inewlam).^2); 
    end
    
    rlamhist2(i) = rlam(id);    
    ilamhist2(i) = ilam(id);  
    
    % -------------        
    if useSBP
        [~, D] = eig(QSBP+S,1i*MSBP/2);
        lam = diag(D) / size(S,2);
        rlam = real(lam);
        ilam = imag(lam);
        
        if i < 3
            [~,id] = min(abs(rlam-ell));
        else % find closest linear extrapolation assuming dL constant
            newlam = rlamhistSBP(i-2) + (rlamhistSBP(i-1) - rlamhistSBP(i-2))*2;
            inewlam = ilamhistSBP(i-2) + (ilamhistSBP(i-1) - ilamhistSBP(i-2))*2;
            [~,id] = min(abs(rlam-newlam).^2 + abs(ilam-inewlam).^2);
        end
        
        rlamhistSBP(i) = rlam(id);
        ilamhistSBP(i) = ilam(id);
    else
        rlamhistSBP(i) = nan;
        ilamhistSBP(i) = nan;
    end
    
end

Lvec = Lvec/pi;
ilamhist = ilamhist/pi;
rlamhist = rlamhist/pi;
ilamhist2 = ilamhist2/pi;
rlamhist2 = rlamhist2/pi;
ilamhistSBP = ilamhistSBP/pi;
rlamhistSBP = rlamhistSBP/pi;

figure(2)
loglog(Lvec,abs(ilamhist),'bo-','linewidth',2,'DisplayName','Full');
hold on
plot(Lvec,abs(ilamhist2),'rs--','linewidth',2,'DisplayName','Lumped');
plot(Lvec,abs(ilamhistSBP),'m^--','linewidth',2,'DisplayName','SBP');
loglog(Lvec,Lvec.^(2*N+2),'b--','DisplayName','h^{2p+2}')
loglog(Lvec,.1*Lvec.^(N+2),'r--','DisplayName','h^{p+2}')
loglog(Lvec,.1*Lvec.^(2*N+1),'k--','DisplayName','h^{2p+1}')

axis on
grid on
set(gca,'fontsize',14)
xlabel('Degrees of freedom per wavelength','fontsize',14)
ylabel('Dissipation relation','fontsize',14)
legend show

figure(3)
loglog(Lvec,abs(rlamhist-Lvec),'bo-','linewidth',2,'DisplayName','Full');
hold on;
plot(Lvec,abs(rlamhist2-Lvec),'rs--','linewidth',2,'DisplayName','Lumped');
plot(Lvec,abs(rlamhistSBP-Lvec),'m^--','linewidth',2,'DisplayName','SBP');
loglog(Lvec,Lvec.^(2*N+2),'b--','DisplayName','h^{2p+2}')
loglog(Lvec,.1*Lvec.^(N+2),'r--','DisplayName','h^{p+2}')
loglog(Lvec,5*Lvec.^(2*N+1),'k--','DisplayName','h^{2p+1}')
axis on
grid on
set(gca,'fontsize',14)
xlabel('Degrees of freedom per wavelength','fontsize',14)
ylabel('Dispersion relation','fontsize',14)

legend show