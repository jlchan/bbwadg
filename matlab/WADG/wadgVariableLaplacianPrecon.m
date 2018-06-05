clear
Globals1D;

f = @(x) (x < -pi/8) + (x > pi/8).*exp(x).*(sin(1+pi*x));
f = @(x) (x > .1) + exp(x).*sin(1+pi*x);

c = 1e4;
k = 1;
w = @(x) c*(1 + .01 + (abs(x)<.5) + sin(1+k*pi*x));
K1D = 1;

N = 7;

r = JacobiGL(0,0,N);

[rq wq] = JacobiGL(0,0,N);

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

StartUp1D;

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

W = w(rq);
plot(r,W,'o--');%return

bmask = [0;ones(size(r,1)-2,1);1];
F = spdiag(bmask);
B = spdiag(1-bmask);


Kw = M + Dr'*Vq'*diag(wq.*W)*Vq*Dr;
K = Dr'*M*Dr;
Kinvw = Dr'*Vq'*diag(wq./w(rq))*Vq*Dr;

Kw = F*Kw*F + B;
K = F*K*F + B;
Kinvw = F*Kinvw*F + B;

invKwadg = K\(Kinvw/K);
W = w(xq); fprintf('contrast = %f\n',max(W)./min(W))

cond(K\Kw)
cond(invKwadg*Kw)
[W D] = eig(inv(invKwadg),Kw); lam = diag(D);
plot(lam,'o')

u = rand(N+1,1);
% hold on
% plot(eig(Kw),'x')



%%

M = 500;

x = linspace(-1 + 2/M,1,M);
[X Y] = meshgrid(x,x);
g = @(x,y)(sin(x*pi).*sin(y*pi));

f  = g( X, Y );


f_f = fft2(f);
u = zeros(size(f_f));
four_coef = [0:M/2 (-M/2+1):-1];

for i = 1:numel(four_coef)
    for j = 1:numel(four_coef)
        if i+j > 2
            u(i,j) = -f_f(i,j)/(pi^2*(four_coef(i)^2+four_coef(j)^2));
        end
    end
end

u = ifft2(u);
figure(5);
mesh(X, Y, real(u));




%% 2D

clear
Globals2D

k = 1;
c = 1e1;
w = @(x,y) c*(1 + 1/c + sin(.5+k*pi*x).*sin(.5+k*pi*y));

K1D = 1;
N = 25;

%[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
[VX VY] = Nodes2D(1); K = 1; EToV = 1:3;


StartUp2D;

Nq = 2*N+1;

[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;

W = w(xq,yq); fprintf('contrast = %f\n',max(W)./min(W))
% W = ones(size(xq));
% ids = randi(length(xq(:)),round(length(xq(:))/3),1);
% W(ids) = 10000;
% color_line3(xq,yq,W,W,'.'); return

bmask = zeros(Np,1); bmask(Fmask(:,1)) = 1;
F = diag(1-bmask); B = diag(bmask);

K = Dr'*M*Dr + Ds'*M*Ds;
Kw = Dr'*(Vq'*diag(wq.*W)*Vq)*Dr + Ds'*(Vq'*diag(wq.*W)*Vq)*Ds;
Kinvw = Dr'*(Vq'*diag(wq./W)*Vq)*Dr + Ds'*(Vq'*diag(wq./W)*Vq)*Ds;

K = F*K*F+B;
Kw = F*Kw*F+B;
Kinvw = F*Kinvw*F+B;

invKwadg = K\(Kinvw/K);

% cond(Kw)
cond(K\Kw)
cond(invKwadg*Kw)

%% larger mesh


clear
Globals2D

k = 1;
c = 1e7;
w = @(x,y) c*(1 + 1/c + sin(.5+k*pi*x).*sin(.5+k*pi*y));

K1D = 4;
N = 10;

% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
[VX VY] = Nodes2D(1); K = 1; EToV = 1:3;

StartUp2D;

Nq = 2*N+1;

[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;

W = w(xq,yq); %fprintf('contrast = %f\n',max(W)./min(W))
W = W(:);

R = getCGRestriction();

VQ = kron(speye(K),Vq);
wJq = diag(wq)*(Vq*J); 
wJq = wJq(:);
Dx = kron(spdiag(rx(1,:)),Dr) + kron(spdiag(sx(1,:)),Ds);
Dy = kron(spdiag(ry(1,:)),Dr) + kron(spdiag(sy(1,:)),Ds);
MM = kron(spdiag(J(1,:)),M);

KK = Dx'*MM*Dx + Dy'*MM*Dy;
KK(vmapB,:) = 0; KK(:,vmapB) = 0; KK(vmapB,vmapB) = eye(length(vmapB)); 
Kh = R*KK*R';

Kw = Dx'*(VQ'*diag(wJq.*W)*VQ)*Dx + Dy'*(VQ'*diag(wJq.*W)*VQ)*Dy;
Kw(vmapB,:) = 0; Kw(:,vmapB) = 0; Kw(vmapB,vmapB) = eye(length(vmapB)); 
Kw = R*Kw*R';

Kinvw = Dx'*(VQ'*diag(wJq./W)*VQ)*Dx + Dy'*(VQ'*diag(wJq./W)*VQ)*Dy;
Kinvw(vmapB,:) = 0; Kinvw(:,vmapB) = 0; Kinvw(vmapB,vmapB) = eye(length(vmapB)); 
Kinvw = R*Kinvw*R';

invKwadg = Kh\(Kinvw/Kh);

% cond(Kw)
cond(Kh\Kw)
cond(invKwadg*Kw)

%% cheb poly

N = 10;
r = JacobiGL(0,0,N);

[rq wq] = JacobiGQ(0,0,N+10);

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

StartUp1D;

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

K = Dr'*M*Dr;
K(1,:) = 0;K(:,1) = 0;K(1,1) = 1;
K(N+1,:) = 0;K(:,N+1) = 0;K(N+1,N+1) = 1;
T = @(j,x) cos(j*acos(x));
for j = 0:N
    V(:,j+1) = T(j,r);
    Vq(:,j+1) = T(j,rq);
end
Dr = V\Dr*V;

% K = Dr'*Vq'*diag(wq)*Vq*Dr;

K = V'*K*V;
