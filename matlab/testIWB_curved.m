function testIWB_curved

clear -global *

Globals3D;

% Order of polymomials used for approximation
N = 2;

filename = 'Grid/sphere48.msh';
filename = 'Grid/sphere400.msh';
% filename = 'Grid/sphere384.msh';
% filename = 'Grid/sphere2736.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

%%

[rp sp tp] = EquiNodes3D(50); %[rp sp tp] = xyztorst(rp,sp,tp);
rf = [r(Fmask(:,1)); r(Np)];
sf = [s(Fmask(:,1)); s(Np)];
tf = [t(Fmask(:,1)); t(Np)];
VWBf = VandermondeGH(N,rp,sp,tp)/VandermondeGH(N,rf,sf,tf);

plot3(rf,sf,tf,'o')

xf = rf ;
yf = sf ;
zf = tf - .25*rand(size(rf));

xyz = VWBf*[xf yf zf];
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

h = color_line3(x,y,z,z,'.');
set(h,'markersize',16)


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




% returns vandermonde matrix and ids of vertex/edge/faces shape functions
function [V v_ids e_ids f_ids] = JVandermonde3D(N, r, s, t)

% vertex ordering
% 3
% |  \
% | 4 \
% 1 -- \ 2
%


Nfp = 4 + 6*(N-1) + 4*max(0,N-2);
V = zeros(length(r), Nfp);

% define linear vertex shape functions
V(:,1) = -(1+r+s+t)/2;
V(:,2) = (1+r)/2;
V(:,3) = (1+s)/2;
V(:,4) = (1+t)/2;
sk = 5;

edges = [1 4; 2 4; 3 4; 1 2; 2 3; 3 1];
e_ids = [];
for e = 1:6            
    for j=0:N-2
        v1 = edges(e,1);     v2 = edges(e,2);
        t = V(:,v1)-V(:,v2);
        V(:,sk) = V(:,v1).*V(:,v2).*JacobiP(t, 1, 1, j);
        e_ids = [e_ids sk];
        sk = sk+1;
    end
end

faces = [1 2 4; 2 3 4; 3 1 4; 1 2 3];
f_ids = [];
for f = 1:4
    v1 = faces(f,1); v2 = faces(f,2); v3 = faces(f,3);
    [x y] = eqbarytoxy(V(:,v1),V(:,v2),V(:,v3)); [rr ss] = xytors(x,y);
    Vf = Vandermonde2D(N-3,rr,ss);    
    for i = 1:size(Vf,2)
        V(:,sk) = V(:,v1).*V(:,v2).*V(:,v3).*Vf(:,i);
        f_ids = [f_ids sk];
        sk = sk + 1;
    end
end

v_ids = 1:4;

return;