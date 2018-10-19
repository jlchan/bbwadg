% wadg weighted convergence (not curved)

clear
Globals2D

k = 1;
f = @(x,y) exp((x+y)).*sin(k*pi*x).*sin(k*pi*y) + 0*(x+y > sin(pi*x));

Kvec = [4 8 16 32 64];
h = 2./Kvec;

sk = 1;
for K1D = Kvec
    
    
    N = 5;
    
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    StartUp2D;
    
    Nq = 2*N+1;
    
    [rq sq wq] = Cubature2D(Nq);
    Vq = Vandermonde2D(N,rq,sq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq));
    xq = Vq*x; yq = Vq*y;
    
    w = 1 + exp(xq+yq).*sin(pi*xq).*sin(pi*yq);
    
    u = Pq*f(xq,yq);
    
    wJq = diag(wq)*(Vq*J);
    
    fwadg1 = Pq*((Vq*u).*w);
    
    MM = 1;
    d = zeros(Np,1);
    skk = 1;
    for i = 0:N
        for j = 0:N-i
            if (i+j <= MM)
                d(skk) = 1;
            end
            skk = skk + 1;
        end
    end
    F = V*(diag(d)/V);
    w2 = Vq*F*Pq*w;
    fwadg2 = Pq*((Vq*u).*w2);
    
    err = Vq*(fwadg1-fwadg2);
    %     err3 = f(xq,yq)-Vq*fwadg;
    
    v = Vq*exp(x+.5*y);
    L2err_diff(sk) = sqrt(abs(sum(sum(wJq.*err.*v))));
    
    sk = sk + 1;
    
end

%
loglog(h,L2err_diff,'k^--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,.1*h.^(MM+1),'b-','linewidth',2)
% loglog(h,1*h.^(N+2),'b-','linewidth',2)
grid on
