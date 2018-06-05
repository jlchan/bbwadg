clear
Globals2D

N = 4;
K1D = 8;
if 0
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    a = .125;
    VX = VX + a*sin(pi*VX).*sin(pi*VY);
    VY = VY + a*sin(pi*VX).*sin(pi*VY);
else
    % one element
    Nv = 3; [VX VY] = EquiNodes2D(1); K = 1; EToV = 1:3;
end

StartUp2D;

% plot(VX,VY,'o');return

[R vmapBT] = getCGRestriction();

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;
wJq = diag(wq)*Jq; wJq = wJq(:);
xq = xq(:); yq = yq(:);

% Vq = bern_basis_tri(N,rq,sq);
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); 

VqK = kron(speye(K),Vq);
MK = kron(spdiag(J(1,:)),M);
Mh = R*MK*R';

Dx = kron(spdiag(rx(1,:)),Dr);
KK = (Dx'*MK*Dx);
Kh = R*KK*R';

fq = ones(size(xq)); 
fq = xq;
fq = sin(pi*xq).*sin(pi*yq);

% L2 projection
b = (R*VqK'*(wJq.*fq));
u = R'*(Mh\b);
L2err = sqrt(sum(wJq.*(VqK*u - fq).^2))

% quasi-interpolant
Rt = spdiag(1./sum(R,2))*R;
% u = R'*Rt*(MK\(VqK'*(wJq.*fq)));
% L2err = sqrt(sum(wJq.*(VqK*u - fq).^2))

% WADG
Pq = Rt*(kron(speye(K),Pq)); % quasi interpolant
% Pq = Rt*(kron(speye(K),inv(M)))*Rt'*(R*VqK'*kron(speye(K),diag(wq)));
% Pq = (R*kron(eye(K),M)*R')\(R*VqK'*kron(speye(K),diag(wq)));
Mp = Pq*spdiag(1./wJq)*Pq'; 
Mp2 = Rt*kron(spdiag(1./J(1,:)),inv(M))*R'; % block diag
u = R'*(Mp*b);
L2err = sqrt(sum(wJq.*(VqK*u - fq).^2))
[cond(Mh),cond(full(Mp*Mh)),cond(full(Mp2*Mh)),cond(full(Mh))]

% a = .00/K;
% [cond(full(Mp*(Mh+a*Kh))),cond(full(Mh+a*Kh))]
% vv = Vp*reshape(u,Np,K);color_line3(xp,yp,vv,vv,'.')


