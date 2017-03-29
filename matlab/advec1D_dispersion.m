function advec1D_dispersion

NB = 4;
Ksub = 8;
K1D = 1;
N = NB+Ksub-1;
smoothKnots = 1;

global alpha
alpha = 0; % 0 = upwinding

ndofs = (N+1);%*K1D;

r = JacobiGL(0,0,N);
rp = linspace(-1,1,500)';
Vp = bsplineVDM(NB,Ksub,rp,smoothKnots); 
[BVDM M Dr] = bsplineVDM(NB,Ksub,r,smoothKnots); % VDM for interp, mass, M\S

S = M*Dr;

e0 = zeros(ndofs,1); e0(1) = 1;
eN = zeros(ndofs,1); eN(ndofs) = 1;

dL = .025;

Lmin = dL;
Lmax = 8*pi;
Lvec = Lmin:dL:Lmax;
for i = 1:length(Lvec)
    L = Lvec(i);
    scale = 1 / K1D; % scale by 1/K1D to incorporate smaller mesh size
    k = 1i*L*scale; % length of domain
    K = S - alpha/2*eN*(eN' - exp(k)*e0') + (2-alpha)/2*e0*(e0' - exp(-k)*eN');      
    [W D] = eig(K,1i*M/2);
    lam = diag(D) / (scale); % normalize by # dofs
    rlam = real(lam);
    ilam = imag(lam);
    
    if i < 3
        [~,id] = min(abs(rlam-L));
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
% plot(Lvec,(rlamhist),'.');hold on;plot(Lvec,Lvec,'k--')
loglog(Lvec,abs(ilamhist),'x');hold on
% loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'o','DisplayName',sprintf('N = %d, Ksub = %d\n',NB,Ksub));hold on
