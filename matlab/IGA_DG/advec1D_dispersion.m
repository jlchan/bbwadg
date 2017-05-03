% function advec1D_dispersion
clear
NB = 3;
Ksub = 32;
K1D = 1;
N = NB+Ksub-1;
smoothKnots = 0;

global alpha
alpha = 0; % 0 = upwinding

ndofs = (N+1)*K1D

r = JacobiGL(0,0,N);
rp = linspace(-1,1,500)';
Vp = bsplineVDM(NB,Ksub,rp,smoothKnots); 
[BVDM M Dr] = bsplineVDM(NB,Ksub,r,smoothKnots); % VDM for interp, mass, M\S

S = M*Dr;

e0 = zeros(N+1,1); e0(1) = 1;
eN = zeros(N+1,1); eN(N+1) = 1;

h = 1/K1D;

Npts = 100;
Lmin = 1/Npts;
Lmax = pi-1/Npts;
Lvec = linspace(Lmin,Lmax,Npts);
for i = 1:length(Lvec)
    ell = Lvec(i);
    k = 1i*ell*(N+1); % length of domain
    A = S - alpha/2*eN*(eN' - exp(k)*e0') + (2-alpha)/2*e0*(e0' - exp(-k)*eN');      
    [W D] = eig(A,1i*M/2);
    lam = diag(D) / (N+1);
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
%     plot(L,rlam,'o');hold on
end

% keyboard

% plot(Lvec,(ilamhist),'o');hold on;plot(Lvec,zeros(size(Lvec)),'k--')
% figure(1)
% plot(Lvec/pi,(rlamhist)/pi,'o-','linewidth',2);hold on;plot(Lvec/pi,Lvec/pi,'k--','linewidth',2)
loglog(Lvec,abs(ilamhist),'x--');hold on
% axis on
% set(gca,'fontsize',14)
% xlabel('Degrees of freedom per wavelength','fontsize',14)
% ylabel('Dispersion relation','fontsize',14)
% legend('Numerical','Exact')
% return
% hold on
% figure(2)
% loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'o','DisplayName',sprintf('N = %d, Ksub = %d\n',NB,Ksub));hold on
