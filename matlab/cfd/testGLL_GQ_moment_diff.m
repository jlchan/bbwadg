clear
clear -global
Globals1D;

N = 3;
Kvec = 2.^(0:5);

u = @(x) exp(sin(2*x));
du = @(x) exp(sin(2*x)).*2.*cos(2*x);
v = @(x) x.^N + 0*cos(x);
iex = integral(@(x) du(x).*v(x),-1,1,'AbsTol',1e-10);

for kk = 1:length(Kvec)
    K1D = Kvec(kk);
    
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    StartUp1D;
    
    [r w] = JacobiGL(0,0,N);
    [rq wq] = JacobiGQ(0,0,N);
    
    Vq = Vandermonde1D(N,rq)/V;
    M = (Vq'*diag(wq)*Vq);
    Pq = M\(Vq'*diag(wq));
    xq = Vq*x;
    Vf = Vandermonde1D(N,[-1;1])/V;
    xf = Vf*x;
    wJq = diag(wq)*(Vq*J);
    wJ = diag(w)*J;
    
    Qr = diag(w)*Dr;
    Qrq = diag(wq)*Vq*(Dr/Vq);
    
    i1 = 0;
    i2 = 0;
    i3 = 0;
    for e = 1:K
        
        % GLL collocation
        i1 = i1 + J(1,e)*rx(1,e)*v(x(:,e))'*Qr*u(x(:,e));
        
        % GQ collocation
        i2 = i2 + J(1,e)*rx(1,e)*v(xq(:,e))'*Qrq*u(xq(:,e));
        
        % weak GQ
        i3 = i3 + sum([-1;1].*u(xf(:,e)).*v(xf(:,e)))-J(1,e)*rx(1,e)*v(xq(:,e))'*Qrq'*u(xq(:,e));
        %         i2 = i2 + J(1,e)*rx(1,e)*(Vq\v(xq(:,e)))'*Qr*(Vq\u(xq(:,e)));
        %         i2 = i2 + J(1,e)*rx(1,e)*(Vq*v(x(:,e)))'*Qrq*(Vq*u(x(:,e)));
    end
    errGLL(kk) = abs(i1 - iex);
    errGQ(kk) = abs(i2 - iex);
    errGQw(kk) = abs(i3 - iex);
end

h = 1./Kvec;
loglog(h,errGLL,'bo--','linewidth',2)
hold on
loglog(h,errGQ,'rx--','linewidth',2)
loglog(h,errGQw,'ks--','linewidth',2)
legend('GLL colloc','GQ colloc','Weak GQ')
% if mod(N,2)==0
%     loglog(h,1e-1*h.^(N),'b--')
% else
%     loglog(h,1e-1*h.^(N+1),'b--')
% end
loglog(h,5e-2*h.^(2*N),'r--')
loglog(h,5e-3*h.^(2*N+2),'k--')

% % test accuracy
% for kk = 1:length(Kvec);
% end
% uh = 1/J(1,e)*(M\(Vf'*diag([-1;1])*u([-1;1])-J(1,e)*rx(1,e)*(Vq*Dr)'*diag(wq)*u(xq(:,e))));
% v = (Pq*du(xq(:,e)) - uh);
% v'*M*v
