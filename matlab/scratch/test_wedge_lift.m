% function test_BB_wedge

N = 3;

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

[re se te] = wedge_nodes(N);
% plot3(re,se,te,'.','markersize',32);return

if 0
    Nplot = 75;
    t1D = linspace(-1,1,Nplot); t1D = t1D(:);
    [rt st] = EquiNodes2D(Nplot); [rt st] = xytors(rt,st); rt = rt(:); st = st(:);
    
    [r t] = meshgrid(rt,t1D);
    [s, ~] = meshgrid(st,t1D);
    
    %V = wedge_sem_basis(N,r,s,t);
    V = bern_wedge(1,r,s,t);
    
    [re se te] = wedge_nodes(1);
    Np = length(re);
    for j = 1:size(V,2)
        text(re-.1,se,te,num2str((1:Np)'))
        color_line3(r,s,t,V(:,j),'.')
        view(3)
        pause
    end
    return
end


%% face quadrature

quadFace = 1;
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
    fids = find(abs(te+1)<1e-8);
end

% plot3(rf,sf,tf,'.','markersize',32); return

%% make geom mapping

[VX VY VZ] = wedge_nodes(1); EToV = 1:6; K = 1;
VZ(1) = VZ(1) + .25;
VZ(6) = VZ(6) + .25;
% VX = VX + .25*randn(size(VX));
% VY = VY + .25*randn(size(VY));
% VZ = VZ + .25*randn(size(VY));

[x y z J geo] = wedge_map(r,s,t,VX,VY,VZ);
[xf yf zf Jf geof] = wedge_map(rf,sf,tf,VX,VY,VZ);
rx = geo.rx; ry = geo.ry; rz = geo.rz;
sx = geo.sx; sy = geo.sy; sz = geo.sz;
tx = geo.tx; ty = geo.ty; tz = geo.tz;

rxf = geof.rx; ryf = geof.ry; rzf = geof.rz;
sxf = geof.sx; syf = geof.sy; szf = geof.sz;
txf = geof.tx; tyf = geof.ty; tzf = geof.tz;

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
if 1
    h = color_line3(r,s,t,J,'.'); set(h,'markersize',32);
%     hold on;
%     h2 = color_line3(rf,sf,tf,sJ,'.'); set(h2,'markersize',64)
%     quiver3(xf,yf,zf,nx,ny,nz)
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
% return


