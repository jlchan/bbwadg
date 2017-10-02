function sphere_gordon_hall_example

clear -global *

Globals3D;

% Order of polymomials used for approximation
N = 3;

filename = 'Grid/sphere48.msh';
% filename = 'Grid/sphere385.msh';
% filename = 'Grid/sphere10970.msh';
% filename = 'Grid/sphere384.msh';
% filename = 'Grid/sphere2736.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

% plot3(r,s,t,'o');text(r+.1,s,t,num2str((0:length(r)-1)'))

%%
if 0
    hold on
    for e = 1:K
        vv = EToV(e,:);
        vx = VX(vv);
        vy = VY(vv);
        vz = VZ(vv);
        %         plot3(vx,vy,vz,'o-');
        if ismember(e,EToE(e,:))
            r = vx.^2 + vy.^2 + vz.^2;
            bids = abs(r-1)<1e-8;
            if (nnz(bids)<3)
                keyboard
            end
            %         plot3(vx(bids),vy(bids),vz(bids),'o','markersize',16);
        end
    end
    
    %     vp = [1 2 3 1 2 4 1 2 3 4 3 1 4 1];
    %     plot3(VX(vv(vp)),VY(vv(vp)),VZ(vv(vp)),'*-');
    keyboard
    
end


%%

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
[Vrq Vsq Vtq] = GradVandermonde3D(N,rq,sq,tq);
Drq = Vrq/V;  Dsq = Vsq/V;  Dtq = Vtq/V;

[rfq sfq tfq wfq] = tet_surface_cubature(30);
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
[Vrfq Vsfq Vtfq] = GradVandermonde3D(N,rfq,sfq,tfq);
Drfq = Vrfq/V;  Dsfq = Vsfq/V;  Dtfq = Vtfq/V;

Nfq = size(rfq,1)/4;

% build W&B operator for a single face
oppvert = [Np (N+1)*(N+2)/2 1 N+1];
rf = [r(Fmask(:,1)); r(Np)];
sf = [s(Fmask(:,1)); s(Np)];
tf = [t(Fmask(:,1)); t(Np)];
VWBf = VandermondeGH(N,r,s,t)/VandermondeGH(N,rf,sf,tf);

VWBe = VandermondeGH(N,r,s,t)/VandermondeGH(N,rf,sf,tf);
% plot3(rf,sf,tf,'o');return

%%

% extract face nodes on the boundary
xf = [];yf = [];zf = [];
xb = [];yb = [];zb = [];
bface = []; eBdry = [];
for e = 1:K
    for f = 1:Nfaces
        if (EToE(e,f)==e)
            xf = [xf x(Fmask(:,f),e)];
            yf = [yf y(Fmask(:,f),e)];
            zf = [zf z(Fmask(:,f),e)];
            xb = [xb x(oppvert(f),e)];
            yb = [yb y(oppvert(f),e)];
            zb = [zb z(oppvert(f),e)];
            bface = [bface f]; % face that's on the boundary
            eBdry = [eBdry e]; % boundary-touching elements
        end
    end
end

% get polar coords
rad = sqrt(xf.^2 + yf.^2 + zf.^2);
theta = atan2(yf,xf);
phi = acos(zf./rad);

% sphere-mapped points
rad = 1;
xfc = rad.*cos(theta).*sin(phi);
yfc = rad.*sin(theta).*sin(phi);
zfc = rad.*cos(phi);

% blend mapped points into the interior
for e = 1:size(xfc,2)
    
    % set face + one boundary node
    xf  = [xfc(:,e); xb(e)];
    yf  = [yfc(:,e); yb(e)];
    zf  = [zfc(:,e); zb(e)];
    xyz = [xf yf zf];
    
    xyzc = VWBf*xyz;
    
    % this assumes boundary face f==1
    x(:,eBdry(e)) = xyzc(:,1);
    y(:,eBdry(e)) = xyzc(:,2);
    z(:,eBdry(e)) = xyzc(:,3);
    
end

% % map global nodes
% [rx,sx,tx,ry,sy,ty,rz,sz,tz,Jq] = tet_geom_factors(x,y,z,Drq,Dsq,Dtq);
% xq = Vq*x; yq = Vq*y; zq = Vq*z;
%
% xfq = Vfq*x; yfq = Vfq*y; zfq = Vfq*z;
% [nx,ny,nz,sJ] = tet_sgeofacs(x,y,z,Drfq,Dsfq,Dtfq);

%% plot mesh

DTf = delaunay(r(Fmask(:,1)),s(Fmask(:,1))); % triangulate face
DTf([8 10],:) = [];
DTf(end+1,:) = [7 8 11];
DTf(end+1,:) = [8 12 11];

hold on

% sk = 1;
% fv = [1 2 3 4; 1 2 4 1; 2 3 4 2; 3 1 4 1];
% FACES = [];
% for e = 1:K
%     for f = 1:4
%         FACES(sk,:) = EToV(e,fv(f,:));
%         sk = sk + 1;
%     end
% end
% VTX = [VX(:) VY(:) VZ(:)];
% p1 = trimesh(FACES,VTX(:,1),VTX(:,2),VTX(:,3));
% set(p1,'EdgeColor','k');
% % set(p1,'LineWidth',2);
% set(p1,'FaceColor','interp');
% material shiny
% set(p1,'FaceLighting','gouraud');

% hold on;plot3(x,y,z,'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
[rfp sfp] = EquiNodes2D(15); [rfp sfp] = xytors(rfp,sfp);
Vfp = Vandermonde2D(N,rfp,sfp)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
DTp = delaunay(rfp,sfp); % triangulate face

% plot(r(Fmask(:,1)),s(Fmask(:,1)),'o');hold on
% text(r(Fmask(:,1)),s(Fmask(:,1)),num2str((1:size(Fmask,1))'))
% trimesh(DTf,r(Fmask(:,1)),s(Fmask(:,1))); return

Emask1 = find(abs(s(Fmask(:,1))+1)<1e-8);
Emask2 = find(abs(r(Fmask(:,1))+1)<1e-8);
Emask3 = find(abs(r(Fmask(:,1))+s(Fmask(:,1)))<1e-8);
Emask = [Emask1(:),Emask2(:),Emask3(:)];
rp = linspace(-1,1,250)';
Vp1D = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N));
for e = 1:K
    for f = 1:4
        if 1 %(e==EToE(e,f))
            
            xf = Vfp*x(Fmask(:,f),e);
            yf = Vfp*y(Fmask(:,f),e);
            zf = Vfp*z(Fmask(:,f),e);
            pt = trimesh(DTp,xf,yf,zf);
            set(pt,'FaceColor','Interp')
            %set(pt,'FaceVertexCData',repmat(C(mod(e,ncolor)+1,:),length(xf),1))
            %         set(pt,'FaceVertexCData',repmat(.5*[1 1 1],length(xf),1))
            %         set(pt,'FaceColor',C(mod(e,5)+1,:));
            
            xf = x(Fmask(:,f),e);
            yf = y(Fmask(:,f),e);
            zf = z(Fmask(:,f),e);
            %         rad = sqrt(xf.^2 + yf.^2 + zf.^2);
            %         theta = atan2(yf,xf);
            %         phi = acos(zf./rad);
            %         xf = 1.01*rad.*cos(theta).*sin(phi);
            %         yf = 1.01*rad.*sin(theta).*sin(phi);
            %         zf = 1.01*rad.*cos(phi);
            %         ptf = trimesh(DTf,xf,yf,zf);
            %         set(ptf,'LineWidth',1)
            %         set(ptf,'EdgeColor','k')
            %         set(ptf,'FaceColor','interp');
            
            for ed = 1:3
                xfp = Vp1D*xf(Emask(:,ed));
                yfp = Vp1D*yf(Emask(:,ed));
                zfp = Vp1D*zf(Emask(:,ed));
                rad = sqrt(xfp.^2 + yfp.^2 + zfp.^2);
                theta = atan2(yfp,xfp);
                phi = acos(zfp./rad);
                xfp = 1.001*rad.*cos(theta).*sin(phi);
                yfp = 1.001*rad.*sin(theta).*sin(phi);
                zfp = 1.001*rad.*cos(phi);
                plot3(xfp,yfp,zfp,'k-','linewidth',3)
            end
        end
    end
end
% shading interp
% material shiny
% camlight

colormap(gray)
view(3)
axis equal;
axis off
axis image
% set(gcf,'Renderer','painters')

% print(gcf,'-dpng','../docs/figs/spherePlanar2.png')
% print(gcf,'-dpng','-r150','../docs/figs/sphereCurved2.png')
return
keyboard

%% check matching of face nodes

for ee = 1:length(eBdry)
    e1 = eBdry(ee);
    for f1 = 1:4
        clf
        
        e2 = EToE(e1,f1);
        f2 = EToF(e1,f1);
        fids1 = (1:Nfq)' + (f1-1)*Nfq;
        fids2 = (1:Nfq)' + (f2-1)*Nfq;
        
        fids1 = 1:4*Nfq;
        fids2 = fids1;
        xf = xfq(fids1,e1);
        yf = yfq(fids1,e1);
        zf = zfq(fids1,e1);
        plot3(xf,yf,zf,'o')
        hold on
        
        xf = xfq(fids2,e2);
        yf = yfq(fids2,e2);
        zf = zfq(fids2,e2);
        plot3(xf,yf,zf,'*')
        
        pause
    end
end

return

%% check normals and reproduction of linears

for e = 1:K
    for f = 1:4
        if (EToE(e,f)==e)
            fids = (1:Nfq)' + (f-1)*Nfq;
            plot3(xfq(fids,e),yfq(fids,e),zfq(fids,e),'.','markersize',16);hold on
            quiver3(xfq(fids,e),yfq(fids,e),zfq(fids,e),nx(fids,e),ny(fids,e),nz(fids,e))
        end
    end
end

u = x+y+z;
dudx = rx.*(Drq*u) + sx.*(Dsq*u) + tx.*(Dtq*u);

return

%% plot shiny pics

[rt st] = EquiNodes2D(50); [rt st] = xytors(rt,st);
Vfp = Vandermonde2D(N,rt,st)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
% Vfp = eye(Nfp);

C = hsv(6);
clf
hold on
for e = 1:size(xfc,2)
    %     plot3(xfc(:,e),yfc(:,e),zfc(:,e),'o')
    xfp = Vfp*xfc(:,e);    yfp = Vfp*yfc(:,e);    zfp = Vfp*zfc(:,e);
    tri = delaunay(xfp,yfp);
    %     colormat = repmat(C(mod(e,4)+1,:),size(xfp,1),1);
    h = trisurf(tri,xfp,yfp,zfp,mod(e,6)+1);
    plot3(xfc(:,e),yfc(:,e),zfc(:,e),'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
    set(h,'edgecolor','none')
end
% shading interp
% shading faceted
% lighting gouraud
material shiny
camlight

axis equal
view(45,20)
axis off
% keyboard
return

% Gordon Hall blending VDM for face = 1
function V = VandermondeGH(N,r,s,t)

V(:,1) = -(1+r+s+t)/2;
V(:,2) = (1+r)/2;
V(:,3) = (1+s)/2;
V(:,4) = (1+t)/2;

sk = 5;
for i = 0:N-2 % edge basis
    V(:,sk) = V(:,1).*V(:,2).*JacobiP(V(:,1)-V(:,2),0,0,i); sk = sk + 1;
    V(:,sk) = V(:,2).*V(:,3).*JacobiP(V(:,2)-V(:,3),0,0,i); sk = sk + 1;
    V(:,sk) = V(:,3).*V(:,1).*JacobiP(V(:,3)-V(:,1),0,0,i); sk = sk + 1;
end

Nb = (N-3);
Vface = Vandermonde2D(Nb,r,s); % lower degree polynomials
for i = 1:(Nb+1)*(Nb+2)/2 % face nodes
    V(:,sk) = V(:,1).*V(:,2).*V(:,3).*Vface(:,i);
    sk = sk + 1;
end



