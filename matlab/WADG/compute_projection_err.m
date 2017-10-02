clear
Globals2D

N = 4;

sk = 1;
for K1D = [1 2 4 8 16 32]
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    
    StartUp2D;
    
    [rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
    Vp = Vandermonde2D(N,rp,sp)/V;
    xp = Vp*x; yp = Vp*y;
    
    [rq sq w] = Cubature2D(3*N); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    xq = Vq*x; yq = Vq*y;
    Jq = Vq*J;
    
    a = 4;
    cfun = @(x,y) (1 + .5*sin(a*pi*x).*sin(a*pi*y) + (y > 0)); % piecewise smooth velocity
%     cfun = @(x,y) 2 + x; % piecewise smooth velocity
    c = (V*V')*Vq'*(diag(w)*cfun(xq,yq)); % L2 projection - can have discontinuities
%     c = cfun(x,y);
    
    cp = Vp*c;
    %     plot3(xp,yp,cp,'.'); keyboard
    
    a = 2;
    u = sin(a*pi*x).*sin(a*pi*y);
    
    disp(sprintf('min continuous c value = %e\n',min(cp(:))))
    if min(cp(:))<1e-12
        disp('negative velocity')
        keyboard
    end
    
    cq = Vq*c;
    Mref = Vq'*diag(w)*Vq;
    
    for e = 1:K
        Minvc{e} = Vq'*diag(w.*1./cq(:,e))*Vq;
        MinvMcM{e} = Mref*((Vq'*diag(w.*cq(:,e))*Vq)\Mref);
    end
    
    matrix_err = zeros(K,1);
    ww = zeros(Np,K);
    for e = 1:K
        %         D = inv(Minvc{e}) - inv(MinvMcM{e});
        %         [W D] = eig(J(1,e)*D'*Mref*D);
        %         [m_err ids] = sort(diag(D),'descend');
        %         matrix_err(e) = m_err(1);
        %         ww(:,e) = W(:,ids(1));
        
        %         diff = Minvc{e}\(Mref*u) - MinvMcM{e}\(Mref*u);
        %         diff = J(1,e)*diag(w)*(Vq*diff).^2;
        %         matrix_err(e) = sum(diff(:));
        
        lam = eig(Mref\MinvMcM{e});
        matrix_err(e) = min(lam);
    end
    
    err(sk) = max(matrix_err);
    lower_bound(sk) = 1/max(cp(:));
    upper_bound(sk) = 1/min(cp(:));
    sk = sk + 1;
end


h = .5.^(1:length(err));
semilogx(h,err,'o-');hold on;
semilogx(h,lower_bound,'k--','linewidth',2);hold on;
semilogx(h,upper_bound,'k--','linewidth',2);hold on;
C = [h(:).^0 log(h(:))]\log(err(:));
title(sprintf('L2 convergence rate = %f, expected rate %f',C(2),N+1))