% wedge experiments
% clear
N = 3;

[rq sq tq w] = wedge_cubature(N);

% plotting points
[a b c] = meshgrid(linspace(-1,1,25));
a = a(:); b = b(:); c = c(:);
[rp sp tp] = wedge_abctorst(a,b,c);

Vq = wedge_basis(N,rq,sq,tq);

% % vertex nodes
r1 = [-1  1 -1 -1  1 -1]'; s1 = [-1 -1  1 -1 -1  1]'; t1 = [-1 -1 -1  1  1  1]';
a = .25;
VX = r1 + a*randn(size(r1)); 
VY = s1 + a*randn(size(s1)); 
VZ = t1 + a*randn(size(t1)); 

% [r s t] = wedge_nodes(N); 
% [V Vr Vs Vt] = wedge_basis(N,r,s,t);
% Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;
% invV = inv(V);
% Vq = Vq*invV;
% [x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = wedge_geom_factors(VX,VY,VZ,r,s,t);

[xq,yq,zq,~,~,~,~,~,~,~,~,~,Jq] = wedge_geom_factors(VX,VY,VZ,rq,sq,tq);
[Vp Vr Vs Vt] = wedge_basis(N,rp,sp,tp);
[xp,yp,zp,~,~,~,~,~,~,~,~,~,Jp] = wedge_geom_factors(VX,VY,VZ,rp,sp,tp);

f = @(x,y,z) x + 2*y + 3*z;
% u = (Vq'*spdiag(w.*Jq(:))*Vq)\(Vq'*(w.*Jq(:).*f(xq,yq,zq)));
% dudx = rx.*(Dr*u) + sx.*(Ds*u) + tx.*(Dt*u)
% dudy = ry.*(Dr*u) + sy.*(Ds*u) + ty.*(Dt*u)
% dudz = rz.*(Dr*u) + sz.*(Ds*u) + tz.*(Dt*u)

M = Vq'*spdiag(w.*Jq(:))*Vq;
M(abs(M)<1e-6) = 0; 
color_line3(xp,yp,zp,Jp,'.')
return
p = symrcm(M); M = M(p,p); 
P = diag(1./sum(M,2));
% P = diag(1./diag(M));

b = (Vq'*(w.*Jq(:).*f(xq,yq,zq)));
lam = eig(P*M);
Lmax = max(lam);Lmin = min(lam);
[x rvec rhist] = cheb_iter(P*M,P*b,P*b,Lmax,Lmin);
semilogy(rvec,'.-')

%% make face integrals
% wedge face vertex ordering
fvW{1} = [1 3 2]; fvW{2} = [4 5 6];
fvW{3} = [1 2 5 4]; fvW{4} = [3 1 4 6]; fvW{5} = [2 3 6 5];
u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
r1 = v(:); s1 = w(:); t1 = u(:); % flipping coordinates for Gmsh
[rf sf tf wf fids] = surface_cubature(N,r1,s1,t1,fvW);
Nfc = length(wf); 
Vfq = wedge_basis(N,rf,sf,tf);

Mf = Vfq'*spdiag(wf)*Vfq;
keyboard
% spy(M)
% b = (N+1)*(N+2)/2;
% M(1:b,1:b)./M((b+1):2*b,(b+1):2*b)
% M(1:b,1:b)\M((2*b+1):3*b,(2*b+1):3*b)

return
[rt st] = Nodes2D(N);[rt st] = xytors(rt,st);
V = Vandermonde2D(N,rt,st);
[~,w1D] = JacobiGQ(0,0,N);

kappa = cond(kron(diag(1./w1D),V*V')*M);
rate= (sqrt(kappa)-1)/(sqrt(kappa)+1);
rate^5
%%
p = symrcm(M); M = M(p,p);

% NfpTri = (N+1)*(N+2)/2;
% i = 1
% M(1:NfpTri,1:NfpTri)./M(i*NfpTri + (1:NfpTri),i*NfpTri + (1:NfpTri))
%%
[W D] = eig(M);
% color_line3(rp,sp,tp,Vp*W(:,3),'.');view(3);colorbar
NfpTri = (N+1)*(N+2)/2;
i = 1
W(1:NfpTri,1:NfpTri)./W(i*NfpTri + (1:NfpTri),i*NfpTri + (1:NfpTri))
