% checks bern quad lift matrix with variable J

N = 3;

[r1D w1D] = JacobiGQ(0,0,N); 
% [r1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,r1D); w1D = sum(inv(V*V'),2);

[r s] = meshgrid(r1D,JacobiGL(0,0,N)); r = r(:); s = s(:);
[wr ws] = meshgrid(w1D); w = wr(:).*ws(:);

% face cubature points: face 1
rf = r1D; sf = -ones(size(r1D)); wf = w1D;

%% randomly warp points + geofacs
r1 = [-1 1 -1 1]'; s1 = [-1 -1 1 1]';  

% ep = .01; VX = r1 + ep*randn(size(r1)); VY = s1 + ep*randn(size(r1));
VX = r1 + [0 0 0 0]'; 
VY = s1 + [0 1 1 0]';

% eval at cubature points
[V1 V1r V1s] = bern_quad(1,r,s);
x = V1*VX; y = V1*VY;
xr = V1r*VX; xs = V1s*VX; 
yr = V1r*VY; ys = V1s*VY;
J = -xs.*yr + xr.*ys;
rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;

% face cubature pts for normals/ sJ
[V1f V1rf V1sf] = bern_quad(1,rf,sf);
xf = V1f*VX; yf = V1f*VY;
xrf = V1rf*VX; xsf = V1sf*VX; 
yrf = V1rf*VY; ysf = V1sf*VY;
Jf = -xsf.*yrf + xrf.*ysf;
rxf = ysf./Jf; sxf =-yrf./Jf; ryf =-xsf./Jf; syf = xrf./Jf;

% face 1
nx = -sxf; ny = -syf;
sJ = sqrt(nx.^2 + ny.^2);
nx = nx./sJ; ny = ny./sJ;
sJ = sJ.*Jf;

plot(VX([1 2 4 3 1]),VY([1 2 4 3 1]),'o-');hold on
vv = J;
color_line3(x,y,vv,vv,'.')
% quiver(xf,yf,nx,ny)
axis equal
return

%% compute operators

% V1D = bern_basis_1D(N,r1D);
[V Vr Vs] = bern_quad(N,r,s);
M = V'*diag(w.*J)*V;
   
Vf = bern_quad(N,rf,sf); 
Mf = Vf'*diag(wf.*sJ)*Vf;

LIFT = M\Mf;
LIFT(abs(LIFT)<1e-6) = 0; 
LIFTf = LIFT(:,1:N+1);

imagesc(LIFTf)

%%
rp1D = linspace(-1,1,50); [rp sp] = meshgrid(rp1D); rp = rp(:); sp = sp(:);
Vp = bern_quad(N,rp,sp);

vv = Vp(:,3);
color_line3(rp,sp,vv,vv,'.')