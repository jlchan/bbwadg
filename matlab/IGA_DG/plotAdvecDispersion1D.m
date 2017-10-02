% function advec1D_dispersion

clear
smoothKnots = 50;
for Ksub = [1 4 8 16]
    NB = 4;
    K1D = 1;
    N = NB+Ksub-1;
    
    global alpha
    alpha = 0; % 0 = upwinding
    
    ndofs = (N+1)*K1D;
    
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
    
    figure(1)
    if Ksub==1
        plot(Lvec/pi,Lvec/pi,'k--','linewidth',2,'DisplayName','Exact');
        hold on;
        plot(Lvec/pi,(rlamhist)/pi,'o-','linewidth',2,'DisplayName','K=1');
    elseif Ksub==4
        plot(Lvec/pi,(rlamhist)/pi,'x-','linewidth',2,'DisplayName','K=4')
    elseif Ksub==8
        plot(Lvec/pi,(rlamhist)/pi,'^-','linewidth',2,'DisplayName','K=8')
    else
        plot(Lvec/pi,(rlamhist)/pi,'s-','linewidth',2,'DisplayName','K=16')
    end
    
    figure(2)
    if Ksub==1
        %loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'o','DisplayName',sprintf('N = %d, Ksub = %d\n',NB,Ksub));hold on
        loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'o-','linewidth',2,'DisplayName','K=1')
        hold on
    elseif Ksub==4
        loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'x-','linewidth',2,'DisplayName','K=4')
    elseif Ksub==8
        loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'^-','linewidth',2,'DisplayName','K=8')
    elseif Ksub==16
        loglog(Lvec,abs(rlamhist(:)-Lvec(:)),'s-','linewidth',2,'DisplayName','K=16')
    end
    
end
figure(1)
axis on; grid on
set(gca,'fontsize',14)
xlabel('Wavelengths per degree of freedom','fontsize',14)
ylabel('Dispersion relation','fontsize',14)
legend(gca,'show','Location','Best')
figure(2)
axis on; grid on
set(gca,'fontsize',14)
xlabel('Wavelengths per degree of freedom','fontsize',14)
ylabel('Dispersion error','fontsize',14)
hold on;
loglog(Lvec,.1*Lvec.^(2*NB+3),'k--','DisplayName','h^{11}')
legend(gca,'show','Location','Best')

if smoothKnots == 0
    figure(1); print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/dispersionUnif.png')
    figure(2); print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/dispersionErrUnif.png')
else
    figure(1); print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/dispersionSmooth.png')
    figure(2); print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/dispersionErrSmooth.png')
end