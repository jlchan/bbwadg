function Wave3D_wedge

clearvars -global

% Driver script for solving the 3D wave equations
hybridgGlobals3D
hybridgGlobalFlags

% make mesh
wedge_mesh; K = size(EToV,1);
p = [4 5 6 1 2 3]; EToV(:,1:6) = EToV(:,p); % permute wedges so jacobian > 0

% ids = abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8 & abs(abs(VZ)-1)>1e-8; % perturb interior verts
% VZ(ids) = VZ(ids) + .1;

a = .2;
VX = VX + a*randn(size(VX));
VY = VY + a*randn(size(VX));
VZ = VZ + a*randn(size(VX));

N = 3;
hybridgStartUp

global rx ry rz sx sy sz tx ty tz J nx ny nz sJ
global mapM mapP mapB
global Dr Ds Dt LIFT Vq Pq Vf

% EToF = EToF(:,1:5); % trim for wedges
% EToE = EToE(:,1:5); % trim for wedges

% drawMesh;pause;close
% plot3(VX,VY,VZ,'o')
% keyboard

%% construct nodal wedges
Nptri = (N+1)*(N+2)/2;
Npquad = (N+1)^2;
Nfp = 2*Nptri + 3*Npquad;

NODETOL=1e-6;
% [r s t] = wedge_nodes(N);

% tensor product wedge nodes
[rtri stri] = Nodes2D(N); [rtri stri] = xytors(rtri,stri);
r1D = JacobiGQ(0,0,N);
[rquad squad] = meshgrid(r1D,r1D); 
rquad = rquad(:); squad = squad(:);
[~,r] = meshgrid(r1D,rtri); [s,t] = meshgrid(r1D,stri);
r = r(:); s = s(:); t = t(:);

keyboard
plot3(r,s,t,'o');return
Np = length(r);

offset = 0;
fids{1} = (1:Nptri) + offset; offset = offset + Nptri;
fids{2} = (1:Nptri) + offset; offset = offset + Nptri;
fids{3} = (1:Npquad) + offset; offset = offset + Npquad;
fids{4} = (1:Npquad) + offset; offset = offset + Npquad;
fids{5} = (1:Npquad) + offset; 

rf = zeros(Nfp,1);
sf = zeros(Nfp,1);
tf = zeros(Nfp,1);

% s = 1
rf(fids{1}) = rtri;
sf(fids{1}) = ones(Nptri,1);
tf(fids{1}) = stri;

% s = -1
rf(fids{2}) = rtri;
sf(fids{2}) = -ones(Nptri,1);
tf(fids{2}) = stri;

% t = -1
rf(fids{3}) = rquad;
sf(fids{3}) = squad;
tf(fids{3}) = -ones(Npquad,1);

% r = -1
rf(fids{4}) = -ones(Npquad,1);
sf(fids{4}) = rquad;
tf(fids{4}) = squad;

% r + t = 0
rf(fids{5}) = rquad;
sf(fids{5}) = squad;
tf(fids{5}) = -rquad;

% plot3(rf,sf,tf,'o');return

x = zeros(Np,K);  y = zeros(Np,K);  z = zeros(Np,K);
rx = zeros(Np,K); sx = zeros(Np,K); tx = zeros(Np,K);
ry = zeros(Np,K); sy = zeros(Np,K); ty = zeros(Np,K);
rz = zeros(Np,K); sz = zeros(Np,K); tz = zeros(Np,K);
J = zeros(Np,K);

nx = zeros(Nfp,K); ny = zeros(Nfp,K); nz = zeros(Nfp,K);
Jf = zeros(Nfp,K);
for e = 1:K
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v); Nv(e) = NvK;
    [xe,ye,ze,rxe,sxe,txe,rye,sye,tye,rze,sze,tze,Je] = wedge_geom_factors(VX(v),VY(v),VZ(v),r,s,t);
    x(:,e) = xe; y(:,e) = ye;  z(:,e) = ze;   
    rx(:,e) = rxe;    sx(:,e) = sxe;    tx(:,e) = txe;
    ry(:,e) = rye;    sy(:,e) = sye;    ty(:,e) = tye;
    rz(:,e) = rze;    sz(:,e) = sze;    tz(:,e) = tze;
    J(:,e) = Je;
        
    [xfe,yfe,zfe,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jfe] = wedge_geom_factors(VX(v),VY(v),VZ(v),rf,sf,tf);
    
    offset = 0;    
    
    % s = 1    
    fid = 1:Nptri;
    nx(fid + offset,e) = sxf(fids{1});
    ny(fid + offset,e) = syf(fids{1});
    nz(fid + offset,e) = szf(fids{1});
    offset = offset + Nptri;
    
    % s = -1
    fid = 1:Nptri;
    nx(fid + offset,e) = -sxf(fids{2});
    ny(fid + offset,e) = -syf(fids{2});
    nz(fid + offset,e) = -szf(fids{2});
    offset = offset + Nptri;    
    
    % t = -1
    fid = 1:Npquad;
    nx(fid + offset,e) = -txf(fids{3});
    ny(fid + offset,e) = -tyf(fids{3});
    nz(fid + offset,e) = -tzf(fids{3});
    offset = offset + Npquad;        
    
    % r = -1
    nx(fid + offset,e) = -rxf(fids{4});
    ny(fid + offset,e) = -ryf(fids{4});
    nz(fid + offset,e) = -rzf(fids{4});
    offset = offset + Npquad;
    
    % r+t = 0
    nx(fid + offset,e) = rxf(fids{5})+txf(fids{5});
    ny(fid + offset,e) = ryf(fids{5})+tyf(fids{5});
    nz(fid + offset,e) = rzf(fids{5})+tzf(fids{5});
    offset = offset + Npquad;
        
    Jf(:,e) = Jfe;    
end
sJ = sqrt(nx.^2 + ny.^2 + nz.^2);    
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*Jf([fids{:}],:);

if min(J(:)) < 1e-12
    disp(sprintf('min J = %f\n',min(J(:))))
    return
end

% interp to surface nodes
Vf = [];
for f = 1:5
    Vf = [Vf ; wedge_basis(N,rf(fids{f}),sf(fids{f}),tf(fids{f}))/wedge_basis(N,r,s,t)];
end

% make surface nodes
xf = Vf*x; 
yf = Vf*y;
zf = Vf*z;

% drawMesh;hold on
% plot3(x,y,z,'o')
% return

% e = 1;
% drawMesh(e)
% plot3(xf(:,e),yf(:,e),zf(:,e),'o')
% hold on
% quiver3(xf(:,e),yf(:,e),zf(:,e),nx(:,e),ny(:,e),nz(:,e))
% return

Nfaces = 5;
nodeids = reshape(1:K*Nfp, Nfp, K);
mapM   = nodeids; % initialize to node ids
mapP   = zeros(Nfp, K);

for e=1:K       
       
    for f = 1:Nfaces        
        enbr = EToE(e,f);
        if (enbr==0) % boundary node
            mapP(nodeids(fids{f},e)) = mapM(nodeids(fids{f},e));
        else
            fnbr = EToF(e,f);
                        
            idM = nodeids(fids{f},e); idM = idM(:);
            idP = nodeids(fids{fnbr},enbr); idP = idP(:);            
            
            tmp = ones(1,length(fids{f}));
            xM = xf(idM)*tmp; yM = yf(idM)*tmp; zM = zf(idM)*tmp;
            xP = xf(idP)*tmp; yP = yf(idP)*tmp; zP = zf(idP)*tmp;
            

            % Compute distance matrix            
            D = (xM -xP').^2 + (yM-yP').^2 + (zM-zP').^2;
            
            [idM, idP] = find(abs(D)<NODETOL);
            mapP(fids{fnbr}(idP),enbr) = mapM(fids{f}(idM),e);
        end
        
    end
end
% mapP = mapP(:); mapM = mapM(:);

mapM = reshape(mapM,Nfp,K);
mapP = reshape(mapP,Nfp,K);

% Create list of boundary nodes
mapB = mapM(mapP==mapM);

mapM = reshape(mapM,Nfp,K);
mapP = reshape(mapP,Nfp,K);

%% lift matrices

[V Vr Vs Vt] = wedge_basis(N,r,s,t);
Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;

% quadrature nodes
[r2 t2 wrt] = Cubature2D(2*N+1);
% [s1D ws] = JacobiGQ(0,0,N);
[s1D ws] = JacobiGL(0,0,N);
[sq rq] = meshgrid(s1D,r2); [~, tq] = meshgrid(s1D,t2);
rq = rq(:); sq = sq(:); tq = tq(:);
[wrt ws] = meshgrid(ws,wrt); wq = wrt(:).*ws(:); wq = 4*wq/sum(wq);

Vq = wedge_basis(N,rq,sq,tq)/V;
invM = V*V';
% invM = inv(Vq'*diag(wq)*Vq);
Pq = invM*Vq'*diag(wq);

% lift matrix
E = zeros(Np,Nfp);
off = 0;

% s = 1
Vtri = Vandermonde2D(N,rf(fids{1}),tf(fids{1}));
Mtri = inv(Vtri*Vtri');
Vfq = wedge_basis(N,rf(fids{1}),sf(fids{1}),tf(fids{1}))/V;
E(:,(1:Nptri)+off) = Vfq'*Mtri;
off = off + Nptri;

% s = -1
Vtri = Vandermonde2D(N,rf(fids{2}),tf(fids{2}));
Mtri = inv(Vtri*Vtri');
Vfq = wedge_basis(N,rf(fids{2}),sf(fids{2}),tf(fids{2}))/V;
E(:,(1:Nptri)+off) = Vfq'*Mtri;
off = off + Nptri;

% quad faces
Vquad = Vandermonde2DQuad(N,rf(fids{3}),sf(fids{3}));
Mquad = inv(Vquad*Vquad');
% Mquad = diag(sum(Mquad,2));
Vfq = wedge_basis(N,rf(fids{3}),sf(fids{3}),tf(fids{3}))/V;
E(:,(1:Npquad)+off) = Vfq'*Mquad;
off = off + Npquad;

Vquad = Vandermonde2DQuad(N,sf(fids{4}),tf(fids{4}));
Mquad = inv(Vquad*Vquad');
% Mquad = diag(sum(Mquad,2));
Vfq = wedge_basis(N,rf(fids{4}),sf(fids{4}),tf(fids{4}))/V;
E(:,(1:Npquad)+off) = Vfq'*Mquad;
off = off + Npquad;

Vquad = Vandermonde2DQuad(N,rf(fids{5}),sf(fids{5}));
Mquad = inv(Vquad*Vquad');
% Mquad = diag(sum(Mquad,2)); 
Vfq = wedge_basis(N,rf(fids{5}),sf(fids{5}),tf(fids{5}))/V;
E(:,(1:Npquad)+off) = Vfq'*Mquad;

LIFT = invM*E;

% keyboard

%% check eigs

if 1
    
    U = zeros(Np*K,4);
    R = zeros(Np*K,4);
    A = zeros(4*Np*K);
    for i = 1:4*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        u = reshape(U(:,2),Np,K);
        v = reshape(U(:,3),Np,K);
        w = reshape(U(:,4),Np,K);
        [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w);
        R(:,1) = rhsp(:);
        R(:,2) = rhsu(:);
        R(:,3) = rhsv(:);
        R(:,4) = rhsw(:);
        
        A(:,i) = R(:);
        
        U(i) = 0;
        if(mod(i,round(4*Np*K/10))==0)
            disp(sprintf('on i = %d out of %d\n',i,4*Np*K))
        end
    end
    
    lam = eig(A);
    plot(lam,'o')
    title(sprintf('Max real part = %g\n',max(real(lam))))
    keyboard
end
%% for plotting

[rp tp sp] = wedge_equi_nodes(20);
Vp = wedge_basis(N,rp,sp,tp)/V;
xp = Vp*x; 
yp = Vp*y; 
zp = Vp*z; 

ids = yp > 0;

%% set initial conditions

f = @(x,y,z) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z);
uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);

p = uexf(x,y,z,0);
u = zeros(Np,K); % velocity
v = zeros(Np,K); 
w = zeros(Np,K); 

%%

resp = zeros(Np,K); resu = zeros(Np,K); resv = zeros(Np,K); resw = zeros(Np,K);
rhsp = zeros(Np,K); rhsu = zeros(Np,K); rhsv = zeros(Np,K); rhsw = zeros(Np,K);

% compute time step size
dt = 2/((N+1)*(N+3)*max(Fscale(:)));

FinalTime = 1;

% outer time step loop
figure
time = 0;
tstep = 0;
while time < FinalTime
    
    if(time+dt>FinalTime), 
        dt = FinalTime-time;  
    end

    % low storage RK    
    for INTRK = 1:5        
        [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w);
                                   
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
        w = w + rk4b(INTRK)*resw;                  
    end

    % Increment time
    time = time+dt; tstep = tstep+1;      
    
    if (mod(tstep,10)==0)        
        clf
        pp = Vp*p;
        drawMesh        
        hold on
        color_line3(xp(ids),yp(ids),zp(ids),pp(ids),'.');        
        caxis([-1 1])
        view(30,30);
        xlabel('x'); ylabel('y')
        title(sprintf('time t = %f, max value of u = %i',time,max(abs(p(:)))))
        colorbar
        axis equal
        drawnow
    end           
end

return
% compute error
uexc = uexf(x,y,z,time);
pc = zeros(NcMax,K);
pc(1:NcH,hexK ) = VH*p(1:NpH,hexK );
pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
pc(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
for e = 1:K
    diff = uexc(:,e) - pc(:,e);
    err2 = err2 + wJ(:,e)'*(diff.^2);
end
err = sqrt(err2);
disp(sprintf('L2 err = %f\n',err))
keyboard


function [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w)

hybridgGlobals3D;

global rx ry rz sx sy sz tx ty tz J nx ny nz sJ
global mapM mapP mapB
global Dr Ds Dt LIFT Vq Pq Vf

% fluxes

pf = Vf*p;
uf = Vf*u;
vf = Vf*v;
wf = Vf*w;

dp = pf(mapP)-pf(mapM);
du = uf(mapP)-uf(mapM);
dv = vf(mapP)-vf(mapM);
dw = wf(mapP)-wf(mapM);

dUn = du.*nx + dv.*ny + dw.*nz;

% apply BCs, todo: try p^+ = - p ^ -, which is energy stable
dp(mapB) = -2*pf(mapB);
dUn(mapB) = 0; 

flux_p = .5*(dp - dUn);
flux_u = .5*(dUn - dp);

% volume derivative operators
pr = Dr*p; ps = Ds*p; pt = Dt*p;
dpdx = rx.*pr + sx.*ps + tx.*pt;
dpdy = ry.*pr + sy.*ps + ty.*pt;
dpdz = rz.*pr + sz.*ps + tz.*pt;

ur = Dr*u; us = Ds*u; ut = Dt*u;
vr = Dr*v; vs = Ds*v; vt = Dt*v;
wr = Dr*w; ws = Ds*w; wt = Dt*w;
dudx = rx.*ur + sx.*us + tx.*ut;
dvdy = ry.*vr + sy.*vs + ty.*vt;
dwdz = rz.*wr + sz.*ws + tz.*wt;

divU = dudx + dvdy + dwdz;
rhsp = -divU.*J + LIFT*(flux_p.*sJ);
rhsu = -dpdx.*J + LIFT*(flux_u.*nx.*sJ);
rhsv = -dpdy.*J + LIFT*(flux_u.*ny.*sJ);
rhsw = -dpdz.*J + LIFT*(flux_u.*nz.*sJ);

Jq = Vq*J;
rhsp = Pq*((Vq*rhsp)./Jq);
rhsu = Pq*((Vq*rhsu)./Jq);
rhsv = Pq*((Vq*rhsv)./Jq);
rhsw = Pq*((Vq*rhsw)./Jq);

return

