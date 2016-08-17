function test_wedge_derivs

N = 2;

% vol cubature nodes
[t1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,t1D); w1D = sum(inv(V*V'),2);
% [t1D w1D] = JacobiGQ(0,0,N);

[rt st wt] = Cubature2D(2*N+1);
[r t] = meshgrid(rt,t1D);
[s, ~] = meshgrid(st,t1D);
[w1 w2] = meshgrid(wt,w1D); w = w1(:).*w2(:);


%% plot and label
% % control pts x SEM points
% [rt st] = EquiNodes2D(N); [rt st] = xytors(rt,st); rt = rt(:); st = st(:);
% [te re] = meshgrid(t1D,rt); [~,se] = meshgrid(t1D,st);
% re = re(:); se = se(:); te = te(:);

[r s t] = wedge_interp_nodes(N);
[V Vr Vs Vt] = wedge_sem_basis(N,r,s,t);
Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
Dt = V\Vt; Dt(abs(Dt)<1e-8) = 0;

Nfptri = (N+1)*(N+2)/2;

[rp sp tp] = wedge_equi_nodes(25);
Vp = wedge_sem_basis(N,rp,sp,tp);

% D1 = V\V1; D1(abs(D1)<1e-8) = 0;
% D2 = V\V2; D2(abs(D2)<1e-8) = 0;
% D3 = V\V3; D3(abs(D3)<1e-8) = 0;
% Dr = Dr(1:Nfptri,1:Nfptri);
% Ds = Ds(1:Nfptri,1:Nfptri);
% D1 = D1(1:Nfptri,1:Nfptri);
% D2 = D2(1:Nfptri,1:Nfptri);
% D3 = D3(1:Nfptri,1:Nfptri);
% keyboard

% plot3(re,se,te,'.','markersize',32);return



%% face quadrature

quadFace = 0;
if quadFace
    % quad face
    [rf tf] = meshgrid(JacobiGQ(0,0,N),t1D); % use full quadrature when possible!
    %     [rf tf] = meshgrid(JacobiGQ(0,0,N)); % full quadrature couples LIFT too much....
    sf = -ones(size(rf));
    [wrf wtf] = meshgrid(w1D); wf = wrf(:).*wtf(:);
    fids = find(abs(se+1)<1e-8);
    
else
    % tri face
    [rf sf wf] = Cubature2D(2*N+1);
    tf = [-ones(size(rf)); ones(size(rf))];
    rf = repmat(rf,2,1);
    sf = repmat(sf,2,1);
    wf = repmat(wf,2,1);
    fids = find(abs(t+1)<1e-8);
end

% plot3(rf,sf,tf,'.','markersize',32);

%% make geom mapping

[rt1 st1] = Nodes2D(1); [rt1 st1] = xytors(rt1,st1);
[t1, r1 ] = meshgrid([-1 1],rt1);
[~,s1] = meshgrid([-1 1],st1);
r1 = r1(:); s1 = s1(:); t1 = t1(:);

ep = 3;
% VX = r1 + ep*randn(size(r1)); VY = s1 + ep*randn(size(r1)); VZ = t1 + ep*randn(size(r1));
bcids = find(abs(t1-1)<1e-8 | abs(t1+1)<1e-8) ;
VX = r1; VY = s1; VZ = t1; VZ(bcids) = t1(bcids) + ep*[.21 .35 .1 .25 1 .2]'; % semi-random
% plot3(VX,VY,VZ,'o'); return

[x y z J geo] = wedge_map(r,s,t,VX,VY,VZ);
[xp yp zp Jp geop] = wedge_map(rp,sp,tp,VX,VY,VZ);
rx = geo.rx; sx = geo.sx; tx = geo.tx;
ty = geo.ty; 
rz = geo.rz; sz = geo.sz; tz = geo.tz;

f = @(x,y,z) x.^N + y.^N - 2*z.^N + x.*y;
dfdx = @(x,y,z) N*x.^(N-1) + y;

u = V\f(x,y,z);
norm(Vp*u - f(xp,yp,zp),'fro')

dudx = rx.*(Dr*u) + sx.*(Ds*u) + tx.*(Dt*u);
norm(Vp*dudx - dfdx(xp,yp,zp),'fro')

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);  Vt = Vandermonde2D(N,r2D,s2D); 
% [r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D); Vt = bern_tri(1,r2D,s2D); 

keyboard

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
Nptri = (N+1)*(N+2)/2;
tx_tri = reshape(tz,Nptri,N+1);
% tx_tri = reshape(J,Nptri,N+1);
txp = tx_tri(:,1);

h = color_line3(r2D,s2D,txp,txp,'.'); set(h,'markersize',32)
% h = color_line3(rp,sp,tp,Jp,'.');set(h,'markersize',32)
% h = color_line3(rp,sp,tp,geop.tx,'.');set(h,'markersize',32)


[xf yf zf Jf geof] = wedge_map(rf,sf,tf,VX,VY,VZ);

% face 1
if quadFace
    nx = -sxf; ny = -syf; nz = -szf;
else
    nx = -txf; ny = -tyf; nz = -tzf;
end
sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*Jf;

% plot3(x,y,z,'o','markersize',32);return
if 0
    h = color_line3(r,s,t,J,'.'); set(h,'markersize',32);
    hold on;
    h2 = color_line3(rf,sf,tf,sJ,'.'); set(h2,'markersize',64)
    return
end

%%

% compute lift matrix
V = wedge_sem_basis(N,r,s,t);
M = V'*diag(w.*J)*V;

Vf = wedge_sem_basis(N,rf,sf,tf);
Mf = Vf'*diag(wf.*sJ)*Vf;

LIFT = M\Mf;
LIFT(abs(LIFT)<1e-6) = 0;
LIFT = LIFT(:,fids);
imagesc(LIFT)

keyboard
return

