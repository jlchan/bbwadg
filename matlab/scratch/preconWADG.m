clear
Globals2D

N = 5;
K1D = 1;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% filename = 'Grid/Other/block2.neu';
filename = 'Grid/Maxwell2D/Maxwell1.neu';
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);


StartUp2D;

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

c = @(x,y) 1 + .5*sin(pi*(x+y));
% plot3(xq,yq,w(xq,yq),'.')

M = Vq'*diag(wq)*Vq;
Dx = kron(diag(rx(1,:)),Dr) + kron(diag(sx(1,:)),Ds);
Dy = kron(diag(ry(1,:)),Dr) + kron(diag(sy(1,:)),Ds);

cq = c(xq,yq);
% cq = repmat(cq(1,:),size(xq,1),1);
for e = 1:K    
    Mw{e} = Vq'*diag(wq.*Jq(:,e).*cq(:,e))*Vq;
    Minvw = Vq'*diag(wq.*Jq(:,e)./cq(:,e))*Vq;
    Pw{e} = M\(Minvw/M);    
end
Mw = blkdiag(Mw{:});
Pw = blkdiag(Pw{:});
% M = kron(diag(J(1,:)),M);

cond(Mw)
cond(kron(diag(J(1,:)),M)\Mw)
cond(Pw*Mw)