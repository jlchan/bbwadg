% wedges
clear
N = 2;

VX = [-1 1 -1 -1 1 -1]'; VY = [-1 -1 -1 1 1 1]'; VZ = [-1 -1 1 -1 -1 1]';
% plot3(VX,VY,VZ,'o'); text(VX,VY,VZ,num2str((1:6)'));return

a = .33;
VX = VX+a*randn(length(VX),1);
VY = VY+a*randn(length(VX),1);
VZ = VZ+a*randn(length(VX),1);

% VY(1) = VY(1)+.5;
% plot3(VX,VY,VZ,'o-');return

% tensor product wedge nodes
[rtri stri] = Nodes2D(N); [rtri stri] = xytors(rtri,stri);
r1D = JacobiGQ(0,0,N);
[r1Dq w1Dq] = JacobiGQ(0,0,N);

[~,r] = meshgrid(r1D,rtri); 
[s,t] = meshgrid(r1D,stri);
r = r(:); s = s(:); t = t(:);

Np = length(r);
Nfp = (N+1)^2*3 + (N+1)*(N+2);
[x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = wedge_geom_factors(VX,VY,VZ,r,s,t);

% [rq sq tq wq] = wedge_cubature(N);
[rqtri sqtri wqtri] = Cubature2D(2*N+1);
[~,rq] = meshgrid(r1Dq,rqtri);
[sq,tq] = meshgrid(r1Dq,sqtri);
[wsq,wrtq] = meshgrid(w1Dq,wqtri);
rq = rq(:); sq = sq(:); tq = tq(:);  wq = wsq(:).*wrtq(:);

[rfq tfq sfq wfq fids] = wedge_surface_cubature(N);
wfq = wfq(:);
% plot3(r,s,t,'*')
% plot3(VX,VY,VZ,'o','markersize',10)
% hold on
% plot3(x,y,z,'*');return

% text(r,s,t,num2str((1:length(r))'))
% plot3(rq,sq,tq,'o')

[V Vr Vs Vt] = wedge_basis(N,r,s,t);

% V1(:,1) = (r+t).*(1-s)/4;
% V1(:,2) = (1+r).*(1-s)/4;
% V1(:,3) = (1+t).*(1-s)/4;
% V1(:,4) = (r+t).*(1+s)/4;
% V1(:,5) = (1+r).*(1+s)/4;
% V1(:,6) = (1+t).*(1+s)/4;

[rp sp tp] = wedge_equi_nodes(25);
Vp = wedge_basis(N,rp,sp,tp)/V;
xp = Vp*x; yp = Vp*y; zp = Vp*z;

[Vq Vrq Vsq Vtq] = wedge_basis(N,rq,sq,tq);
Vq = Vq/V; Vrq = Vrq/V; Vsq = Vsq/V; Vtq = Vtq/V;

Vfq = wedge_basis(N,rfq,sfq,tfq)/V;

xf = Vfq*x;
yf = Vfq*y;
zf = Vfq*z;

% vol/surf geofacs
[xq,yq,zq,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq] = wedge_geom_factors(VX,VY,VZ,rq,sq,tq);
[xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = wedge_geom_factors(VX,VY,VZ,rfq,sfq,tfq);
nx = zeros(size(xfq)); ny = zeros(size(xfq)); nz = zeros(size(xfq));
nx(fids{1}) = -sxf(fids{1}); ny(fids{1}) = -syf(fids{1}); nz(fids{1}) = -szf(fids{1});
nx(fids{2}) = -txf(fids{2}); ny(fids{2}) = -tyf(fids{2}); nz(fids{2}) = -tzf(fids{2});
nx(fids{3}) = rxf(fids{3})+txf(fids{3}); 
ny(fids{3}) = ryf(fids{3})+tyf(fids{3}); 
nz(fids{3}) = rzf(fids{3})+tzf(fids{3}); 
nx(fids{4}) = -rxf(fids{4}); ny(fids{4}) = -ryf(fids{4}); nz(fids{4}) = -rzf(fids{4});
nx(fids{5}) = sxf(fids{5}); ny(fids{5}) = syf(fids{5}); nz(fids{5}) = szf(fids{5});
sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*Jf;

nxJ = nx.*sJ;
nyJ = ny.*sJ;
nzJ = nz.*sJ;

% plot3(xfq,yfq,zfq,'o')
% hold on
% quiver3(xfq,yfq,zfq,nx,ny,nz)

Dr = Vr/V; Ds = Vs/V; Dt = Vt/V; 
% f = @(x,y,z) x.^N + y.^N + z.^N;
% dfdx = @(x,y,z) N*x.^(N-1);
% dfdy = @(x,y,z) N*y.^(N-1);
% dfdz = @(x,y,z) N*z.^(N-1);
% norm(Vp*(rx.*(Dr*f(x,y,z)) + sx.*(Ds*f(x,y,z)) + tx.*(Dt*f(x,y,z)))-dfdx(xp,yp,zp),'fro')
% norm(Vp*(ry.*(Dr*f(x,y,z)) + sy.*(Ds*f(x,y,z)) + ty.*(Dt*f(x,y,z)))-dfdy(xp,yp,zp),'fro')
% norm(Vp*(rz.*(Dr*f(x,y,z)) + sz.*(Ds*f(x,y,z)) + tz.*(Dt*f(x,y,z)))-dfdz(xp,yp,zp),'fro')

% check IBP
f = randn(Np,1);

Mref = Vq'*diag(wq)*Vq;
LHSnodal = Mref*(J.*(rx.*(Dr*f) + sx.*(Ds*f) + tx.*(Dt*f)));

% Mf = Vfq'*diag(wfq)*Vfq;
% LHSnodalIBP = Vfq'*diag(wfq.*nx)*Vfq*f - (diag(J.*rx)*Dr' + diag(J.*sx)*Ds' + diag(J.*tx)*Dt')*Mref*f;

LHS = Vq'*(Jq.*wq.*(rxq.*(Vrq*f) + sxq.*(Vsq*f) + txq.*(Vtq*f)));

LHSIBP = Vfq'*(wfq.*nx.*sJ.*(Vfq*f)) - (Vrq'*(rxq.*Jq.*wq.*(Vq*f)) + Vsq'*(sxq.*Jq.*wq.*(Vq*f)) + Vtq'*(txq.*Jq.*wq.*(Vq*f)));

disp(sprintf('diff b/w nodal collocation and quadrature RHS: %g\n',norm(LHSnodal-LHS)))
disp(sprintf('diff b/w quadrature and IBP RHS: %g\n',norm(LHS-LHSIBP)))

rxJ = reshape(rx.*J,(N+1)*(N+2)/2,N+1);
sxJ = reshape(sx.*J,(N+1)*(N+2)/2,N+1);
txJ = reshape(tx.*J,(N+1)*(N+2)/2,N+1);

ryJ = reshape(ry.*J,(N+1)*(N+2)/2,N+1);
syJ = reshape(sy.*J,(N+1)*(N+2)/2,N+1);
tyJ = reshape(ty.*J,(N+1)*(N+2)/2,N+1);

rzJ = reshape(rz.*J,(N+1)*(N+2)/2,N+1);
szJ = reshape(sz.*J,(N+1)*(N+2)/2,N+1);
tzJ = reshape(tz.*J,(N+1)*(N+2)/2,N+1);

J = reshape(J,(N+1)*(N+2)/2,N+1);

r = reshape(r,(N+1)*(N+2)/2,N+1);
s = reshape(s,(N+1)*(N+2)/2,N+1);
t = reshape(t,(N+1)*(N+2)/2,N+1);

% plot3(rfq(fids{2}),sfq(fids{2}),nx(fids{2}).*sJ(fids{2}),'s')
% plot3(rfq(fids{3}),sfq(fids{3}),nx(fids{3}).*sJ(fids{3}),'s')
% plot3(sfq(fids{4}),tfq(fids{4}),nx(fids{4}).*sJ(fids{4}),'s')

Vq = wedge_basis(N,rq,sq,tq)/V;
M = Vq'*diag(wq)*Vq;

% triangular face s = -1
rf = rtri; sf = -ones((N+1)*(N+2)/2,1); tf = stri;
Vft = wedge_basis(N,rf,sf,tf)/V;
Vtri = Vandermonde2D(N,rtri,stri);
Mt = inv(Vtri*Vtri');
Ltri = M\(Vft'*Mt);
Ltri(abs(Ltri)<1e-8) = 0;

% quadrilateral face t = -1
[rquad squad] = meshgrid(r1Dq,r1D); rquad = rquad(:); squad = squad(:);
% plot(rquad,squad,'o');text(rquad(:),squad(:),num2str((1:length(rquad(:)))'));return
rf = rquad; sf = squad; tf = -ones(size(rquad));
Vfq = wedge_basis(N,rf,sf,tf)/V;
Vquad = Vandermonde2DQuad(N,rquad,squad); 
Mq = inv(Vquad*Vquad'); % exact mass matrices

Lquad = M\(Vfq'*Mq);
Lquad(abs(Lquad)<1e-8) = 0;

return
% check metric identities
norm(Dr*rxJ(:) + Ds*sxJ(:) + Dt*txJ(:))
norm(Dr*ryJ(:) + Ds*syJ(:) + Dt*tyJ(:))
norm(Dr*rzJ(:) + Ds*szJ(:) + Dt*tzJ(:))

return
%% pyramids

N = 6;

[VX VY VZ] = pyr_nodes(1);
a = .25;
VX = VX+a*randn(length(VX),1);
VY = VY+a*randn(length(VX),1);
VZ = VZ+a*randn(length(VX),1);
% plot3(VX,VY,VZ,'o')

[r s t] = pyr_nodes(N);
V = pyr_basis(N,r,s,t);

[rq sq tq wq] = pyr_cubature(N);
Vq = pyr_basis(N,rq,sq,tq);
[x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = pyr_geom_factors(VX,VY,VZ,rq,sq,tq);

M = Vq'*diag(wq.*J)*Vq;

[x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = pyr_geom_factors(VX,VY,VZ,r,s,t);
% color_line3(r,s,t,J,'.')

fids = find(abs(s+1)<1e-8);
Vft = pyr_basis(N,r(fids),s(fids),t(fids));
Vtri = Vandermonde2D(N,r(fids),t(fids));
Mt = inv(Vtri*Vtri');

Ltri = M\(Vft'*Mt);
Ltri(abs(Ltri)<1e-8) = 0;

fids = find(abs(t+1)<1e-8);
Vfq = pyr_basis(N,r(fids),s(fids),t(fids));
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
M1D = inv(V1D*V1D'); 
%M1D = diag(sum(M1D,2));
Mq = kron(M1D,M1D);

Lquad = M\(Vfq'*Mq);
Lquad(abs(Lquad)<1e-8) = 0;

% imagesc(abs(invV./min(invV(abs(invV(:))>1e-8))))


