% hermite 1D

% make constraints
Nq = 2; 
N = 2*Nq+2;%(2*Nq-3)*3 + 4;
rc = JacobiGL(0,0,Nq);
V = Vandermonde1D(N,rc);
Vr = GradVandermonde1D(N,rc);

% helpful operators
r = JacobiGL(0,0,N); 
Dr =  Vandermonde1D(N,r)\GradVandermonde1D(N,r); Dr(abs(Dr)<1e-8) = 0;
Vrr = Vandermonde1D(N,rc(2:end-1))*Dr^2; % second derivative

VDM = [V;Vr;Vrr];
invV = inv(VDM);

% make operators
r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r)*invV;
Dr = GradVandermonde1D(N,r)*invV;

rp = linspace(-1,1,250); rp = rp(:);
Vp = Vandermonde1D(N,rp);

ids = 1:N+1; %(1:Nq+1) + Nq+1;
plot(rp,Vp*invV(:,ids))
hold on;plot(rc,0*rc,'o','markersize',24)

[rq w] = JacobiGQ(0,0,N);
Vq = Vandermonde1D(N,rq)*invV;

M = Vq'*diag(w)*Vq;
