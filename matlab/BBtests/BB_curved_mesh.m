Globals2D;

N = 3;

Nq = 2*N+1;

filename = 'Grid/Other/circA01.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% K1D = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% build tri for connectivity
[vx vy] = EquiNodes2D(N); vx = vx(:)'; vy = vy(:)';
tri = delaunayFixArea(vx,vy);
triplot(tri,vx,vy);
h = 2/N; % edge length

N2 = 2*N;
[r2 s2] = Nodes2D(N2); [r2 s2] = xytors(r2,s2);
E = bern_basis_tri(N2,r2,s2)\bern_basis_tri(N,r2,s2);
[vx2 vy2] = EquiNodes2D(N2); vx2 = vx2(:)'; vy2 = vy2(:)';
tri2 = delaunayFixArea(vx2,vy2);

% for mesh
[VX VY] = EquiNodes2D(1); VX = VX(:)'; VY = VY(:)'; K = 1; EToV = 1:3;
StartUp2D;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
[rp sp] = EquiNodes2D(200); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;

x1 = x; y1 = y;
xe1 = Ve*x; ye1 = Ve*y;

opt = 1;
displacement = .25*cos(pi/2*x(Fmask(:,1),:));


% xeC = xe1; yeC = ye1; yeC(3) = yeC(3) + .6;
% displacement = bern_basis_1D(N,JacobiGL(0,0,N))*(yeC(Fmask(:,1)) - ye1(Fmask(:,1)));

% displacement = -.5*(1-x(Fmask(:,1),:)).*(1+x(Fmask(:,1),:));
% displacement = .5*cos(pi/2*x(Fmask(:,1),:));

% y(Fmask(:,1),:) = y(Fmask(:,1),:) + .1*sin(pi*x(Fmask(:,1),:));
% [rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);J = Vp*J;
% [max(J(:)),min(J(:)),max(J(:))/min(J(:))]
% xe1 = Ve*x; ye1 = Ve*y;

[VB VBr VBs] = bern_basis_tri(N,r,s);
DBr = VB\VBr; DBr(abs(DBr)<1e-8) = 0;
DBs = VB\VBs; DBs(abs(DBs)<1e-8) = 0;

xe = Ve*x; ye = Ve*y;
xB = VB\x; yB = VB\y;

% reset interior ids
nfids = setdiff(1:Np,Fmask(:,1));
xB(nfids,:) = xe1(nfids,:);
yB(nfids,:) = ye1(nfids,:);

% totally interior ids
% iids = setdiff(1:Np,Fmask(:)); % anchor all face nodes
iids = nfids; % anchor only bottom face
% iids = setdiff(iids,Np); % anchor top vertex

% spring displacement
if opt==1
    y0 = y(Fmask(:,1),:);
    yF = y0 + displacement;
    yBf = bern_basis_1D(N,JacobiGL(0,0,N))\yF;
    
    % adjacency matrix
    tri = delaunayFixArea(xB',yB');
    A = sparse(Np,Np);
    for e = 1:size(tri,1)
        A(tri(e,:),tri(e,:)) = 1; % adjacency matrix
    end
    %     A = DBr*DBr' + DBs*DBs';    
    A = E'*E;
%         A = double(abs(A)>.1);
    
% keyboard
    nbrs = {};
    weights = {};
    for i = 1:Np
        nbrs{i} = setdiff(find(A(i,:)),i);
        blens = [];
        for jj = 1:length(nbrs{i})
            j = nbrs{i}(jj);
            v = [xe1(i)-xe1(j); ye1(i)-ye1(j)];
            d = norm(v);
            blens = [blens d];
            wvec = A(i,nbrs{i});
            wvec = wvec/max(wvec);
        end
        barlen{i} = blens;
        weights{i} = wvec;
    end
    
    maxit = 200;
    for iter = 1:maxit
        
        alpha = min(1,2*iter/maxit); % give relaxation time
        yB(Fmask(:,1),:) = y0(Fmask(:,1))*(1-alpha) + alpha*yBf;
        
        clf;
        hold on;
%         triplot(tri,xB,yB);plot(xB,yB,'s')
        triplot(tri2,E*xB,E*yB);plot(E*xB,E*yB,'s')
        title(sprintf('iteration %d',iter))
        drawnow
        
        %         [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
        %         Je = Ve*J;
        
        f = zeros(Np,2);
        for i = 1:Np
            
            nbrlist = nbrs{i};
            
            for j = 1:length(nbrlist)
                nbr = nbrlist(j);
                v = [xB(i)-xB(nbr); yB(i)-yB(nbr)];
                d = norm(v);
                h = barlen{i}(j);
%                 k = weights{i}(j);
                k = 1;
                k =  k*min(10,1./abs(h-d))*(d < h);
%                 f(i,:) = f(i,:) +  k*v(:)'*(h-d);
                f(i,:) = f(i,:) + k*v(:)'*(h > d)*(h-d) + .1*v(:)'*(h < d)*(h-d);
            end
            
            
        end
        xB(iids) = xB(iids) + .25*(N/N2)*h*f(iids,1);
        yB(iids) = yB(iids) + .25*(N/N2)*h*f(iids,2);
        
    end
    
    x = VB*xB; y = VB*yB;
end

% W&B
if opt==2
    f = 1;
    %     vr = re; % tangential coord
    vr = r; rv = r; sv = s;% tangential coord
    %     vr = re; rv = re; sv = se; % BB version
    fr = vr(Fmask(:,f));
    
    % evaluate deformation of coordinates
    fdy = displacement;
    %     fdy = bern_basis_1D(N,JacobiGL(0,0,N))\fdy;
    fdx = 0*fdy;
    
    % build 1D Vandermonde matrix for face nodes and volume nodes
    Vface = Vandermonde1D(N, fr);  Vvol  = Vandermonde1D(N, vr);
    % compute unblended volume deformations
    vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);
    
    % blend deformation and increment node coordinates
    ids = find(abs(1-vr)>1e-7); % warp and blend
    if(f==1) blend = -(rv(ids)+sv(ids))./(1-vr(ids)); end;
    if(f==2) blend =      +(rv(ids)+1)./(1-vr(ids)); end;
    if(f==3) blend = -(rv(ids)+sv(ids))./(1-vr(ids)); end;
    %xB(ids) = xB(ids)+blend.*vdx(ids);
    %yB(ids) = yB(ids)+blend.*vdy(ids);
    x(ids) = x(ids)+blend.*vdx(ids);
    y(ids) = y(ids)+blend.*vdy(ids);
    
    %     x = VB*xB; y = VB*yB;
    xB = VB\x; yB = VB\y;
end

% mesh smoothing
if opt==3
    y(Fmask(:,1),:) = y(Fmask(:,1),:) + displacement;
    [rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);
    %     J = Vp*J;
    %     [max(J(:)),min(J(:)),max(J(:))/min(J(:))]
    VB = bern_basis_tri(N,r,s);
    xe = Ve*x; ye = Ve*y;
    xB = VB\x; yB = VB\y;
        
    % simple averaging smoother
    A = zeros(Np);
    for e = 1:size(tri,1)
        A(tri(e,:),tri(e,:)) = 1; % adjacency matrix
    end
    A = abs(A)>1e-8;
    A = A - diag(diag(A)); % remove self interactions
    A = diag(1./sum(A,2))*A;
    
    
    A = VB;   
    
    if 1
        bcids = Fmask(:);
        xBdry = xB(bcids);
        yBdry = yB(bcids);
        for i = 1:50
            xB = A*(xB-xe) + xe;
            yB = A*(yB-ye) + ye;
            xB(bcids) = xBdry;
            yB(bcids) = yBdry;
        end
    end
    %     plot3(xe1,ye1,yB-ye1,'o');return
    x = VB*xB; y = VB*yB;
    
end

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
Jp = Vp*J;

dJ = max(max(abs(Vp*Dr*J)),max(abs(Vp*Ds*J)));
disp(sprintf('opt = %d\n',opt))
[max(Jp(:)),min(Jp(:)),max(Jp(:))/min(Jp(:)),dJ]

clf
%h = color_line3(Vp*x,Vp*y,Jp,Jp,'.'); set(h,'markersize',12)
PlotField2D(50,x,y,J-1.05*max(J(:)));%(xp,yp,'.','markersize',24)
hold on
% plot(x,y,'go','markersize',12,'linewidth',2)

% for i = 1:size(E,2)
%     ids = find(abs(E(:,i))>1e-8);
%     vals = E(ids,i);
%     clf
%     plot(vx2,vy2,'o')
%     hold on
%     plot3(vx2(ids),vy2(ids),vals,'x')
%     view(3)
%     pause
% end

xB2 = E*xB;
yB2 = E*yB;
triplot(tri2,xB2,yB2,'k','linewidth',2);
% hull = convhull(xB2,yB2);
% plot(xB2(hull,:),yB2(hull,:), 'k-','markersize',12,'linewidth',2)
plot3(xB2,yB2,zeros(size(xB2))+.025,'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
view(2)
axis off
axis equal
%  print(gcf,'-dpng','bb_submesh3.png')