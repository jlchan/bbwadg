function sphere_gordon_hall

clear -global *

Globals3D;

% Order of polymomials used for approximation
N = 3;

% filename = 'Grid/sphere48.msh';
filename = 'Grid/sphere385.msh';
% filename = 'Grid/sphere1780.msh';
% filename = 'Grid/sphere10970.msh';
% filename = 'Grid/sphere384.msh';
% filename = 'Grid/sphere2736.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

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

% plot3(x,y,z,'o');return

%% plot mesh

hold on

% hold on;plot3(x,y,z,'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
[rfp sfp] = EquiNodes2D(25); [rfp sfp] = xytors(rfp,sfp);
Vfp = Vandermonde2D(N,rfp,sfp)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
DTp = delaunay(rfp,sfp); % triangulate face

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
%             pt = trimesh(DTp,xf,yf,zf);
            pt = trisurf(DTp,xf,yf,zf);
            set(pt,'FaceColor','Interp')
            set(pt,'edgecolor','none')
            
%             set(pt,'FaceVertexCData',repmat(C(mod(e,ncolor)+1,:),length(xf),1))
%             set(pt,'FaceVertexCData',repmat(.9*[1 1 1],length(xf),1))
            %         set(pt,'FaceColor',C(mod(e,5)+1,:));
            
            xf = x(Fmask(:,f),e);
            yf = y(Fmask(:,f),e);
            zf = z(Fmask(:,f),e);            
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
                plot3(xfp,yfp,zfp,'k-','linewidth',2)
            end
        end
    end
end
shading interp
% material shiny
% camlight

colormap(gray)
view(3)
axis equal;
axis off
axis image
% set(gcf,'Renderer','painters')

% print(gcf,'-dpng','-r150','../spherePlanar1.png')
print(gcf,'-dpng','-r150','../sphereCurved1.png')
return

%% plot shiny pics

[rt st] = EquiNodes2D(25); [rt st] = xytors(rt,st);
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
%     plot3(xfc(:,e),yfc(:,e),zfc(:,e),'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
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



