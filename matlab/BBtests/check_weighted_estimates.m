clear
Globals2D

N = 3;

sk = 1;
for K1D = [4 8 16 32]
    if K1D==1
        [VX VY] = Nodes2D(1); K = 1; EToV = 1:3;
    else
        [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    end
    
    StartUp2D;
    
    [rq sq w] = Cubature2D(3*N); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    xq = Vq*x; yq = Vq*y;
    
    a = 1;
    cfun = @(x,y) (1 + .5*sin(a*pi*x).*sin(a*pi*y) + 1*(y > 0)); % piecewise smooth velocity
%     cfun = @(x,y) 2 + x; % piecewise smooth velocity
    %     c = (V*V')*Vq'*(diag(w)*cfun(xq,yq)); % L2 projection - can have discontinuities
    %     c = cfun(x,y);   % interpolation
    %     cq = Vq*c;
    
    cq = cfun(xq,yq);
    
    Mref = Vq'*diag(w)*Vq;
    
    a = 1;
    f = @(x,y) sin(a*pi*x).*sin(a*pi*y);
    
    wJ = diag(w)*(Vq*J);
    uproj = (V*V')*Vq'*(diag(w)*f(xq,yq));
    
    % solve (1/c * u,v) = (p,v) for u ~ p*c
    for e = 1:K
        M = Vq'*diag(w.*1./cq(:,e))*Vq;
        uweight = M\(Mref*uproj(:,e));
        %diff(:,e) = (Vq*uproj(:,e)).*cq(:,e) -(Vq*(uweight));
        diff(:,e) = f(xq(:,e),yq(:,e)).*cq(:,e) -(Vq*(uweight));
        uprojc(:,e) = uweight;
    end
    diff = wJ.*diff.^2;
    proj_err = sqrt(sum(diff(:)));
    
    % solve (u,v) = (c*Pi*p,v) for u ~ p*c
    for e = 1:K
        M = Vq'*diag(w.*cq(:,e))*Vq;
        Mc = Mref*(M\Mref);
        uweight = Mc\(Mref*uproj(:,e));
        %diff(:,e) = (Vq*uproj(:,e)).*cq(:,e) -(Vq*(uweight));
        %diff(:,e) = f(xq(:,e),yq(:,e)).*cq(:,e) -(Vq*(uweight));
        diff(:,e) = Vq*uprojc(:,e) -Vq*uweight;
    end
    diff = wJ.*diff.^2;
    wproj_err = sqrt(sum(diff(:)));
    
%     % u - c* P_c * u
%     for e = 1:K
%         uweight = (Vq'*diag(w.*cq(:,e))*Vq)\(Mref*uproj(:,e));
%         diff(:,e) = f(xq(:,e),yq(:,e)) - cq(:,e).*(Vq*(uweight));
%     end
%     diff = wJ.*diff.^2;
%     wproj_err = sqrt(sum(diff(:)));            
    
    % u/c - P_c * u
    for e = 1:K
        M = Mref*((Vq'*diag(w.*cq(:,e))*Vq)\Mref);
        uweight = M\(Vq'*(w./cq(:,e).*f(xq(:,e),yq(:,e))));
        diff(:,e) = f(xq(:,e),yq(:,e)) - Vq*(uweight);
    end
    diff = wJ.*diff.^2;
    invwproj_err = sqrt(sum(diff(:)));
    
    err(sk) = proj_err;
    errw(sk) = wproj_err;
    errinvw(sk) = invwproj_err;
    sk = sk + 1;
end

h = 2*.5.^(1:length(err));
loglog(h,err,'o-');
hold on
loglog(h,errinvw,'x-')
loglog(h,errw,'s-') % u - c*P_c*u depends on regularity of c
%legend('1/c-weighted mass error','u/c-P_c(u)','u - cP_c u')
legend('1/c-weighted mass error','u/c-P_c(u)','diff between c-mass and weighted')
% legend('1/c-weighted mass error','c-weighted RHS err')
fit = [log(h(:)) ones(size(h(:)))]\log(err(:));
fitw = [log(h(:)) ones(size(h(:)))]\log(errw(:));
fitinvw = [log(h(:)) ones(size(h(:)))]\log(errinvw(:));

title(sprintf('Order %d: rate = %f, %f, %f\n',N,fit(1),fitw(1) ,fitinvw(1) ))

