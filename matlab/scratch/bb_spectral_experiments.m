%% 2d inverse mass matrix

clear
N = 3;
[r s w] = Cubature2D(2*N);
Vq = bern_basis_tri(N,r,s);
M = Vq'*diag(w)*Vq;
invM = inv(M);
invM = round(invM./min(abs(invM(invM(:)>0))));

r1D = JacobiGL(0,0,N);
for i = 0:N
    for j = i:N
        E{i+1,j+1} = bern_basis_1D(N-i,r1D)\bern_basis_1D(N-j,r1D);
    end
end

%%
clear
N = 2;
[r s t w] = tet_cubature(2*N);
% plot3(r,s,t,'o');
Vq = bern_basis_tet(N,r,s,t);
M = Vq'*diag(w)*Vq;
M = round(M./min(M(:)));

r1D = JacobiGL(0,0,N);
for i = 0:N
    for j = i:N
        E{i+1,j+1} = bern_basis_1D(N-i,r1D)\bern_basis_1D(N-j,r1D);
    end
end

invM = inv(M);
invM = round(invM./min(abs(invM(:))));

%% test ainsworth duffy BB

[ap bp] = meshgrid(linspace(-1,1,50));
ap = ap(:); bp = bp(:);
rp = 0.5*(1+ap).*(1-bp)-1;
sp = bp;
% plot(rp,sp,'o')
% return
VB1 = bern_basis_1D(N,bp);
sk = 1;
for i = 0:N
    VB2 = bern_basis_1D(N-i,ap);
    
    for j = 0:N-i
        clf;
        vv = VB1(:,i+1).*VB2(:,j+1);
        color_line3(rp,sp,vv,vv,'.');
        pause
        sk = sk + 1;
    end
end


%% test lsq fitting of degree reduction
N = 4;
[r w] = JacobiGQ(0,0,N);
V1 = bern_basis_1D(N,r);
V2 = bern_basis_1D(N-1,r);
V3 = bern_basis_1D(N-2,r);

norm(pinv(V1\V3) - pinv(V2\V3)*pinv(V1\V2),'fro') % can do projection sequentially

%% test spectral multiplication

N = 4;
r = linspace(-1,1,N+1)';%JacobiGL(0,0,N);
V = bern_basis_1D(N,r);

f = @(r) exp(r);
g = @(r) sin(r);

hold on
rp = linspace(-1,1,50);
Vp = Vandermonde1D(N,rp)/Vandermonde1D(N,r);
plot(rp,Vp*diag(g(r))*f(r),'^')

Vp = bern_basis_1D(N,rp);
fB = V\f(r); % BB version
T = inv(V)*diag(g(r))*V;
plot(rp,Vp*T*fB,'o')

plot(rp,g(rp).*f(rp))



