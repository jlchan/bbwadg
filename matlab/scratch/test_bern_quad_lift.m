% checks bern quad lift matrix with variable J
function test_bern_quad_lift

N = 3;

% [r1D w1D] = JacobiGQ(0,0,N);
[r1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,r1D); w1D = sum(inv(V*V'),2);

[s r] = meshgrid(r1D); r = r(:); s = s(:);
[wr ws] = meshgrid(w1D); w = wr(:).*ws(:);

% face cubature points:
% face 1
rf = r1D; sf = -ones(size(r1D)); wf = w1D;
% face 2: vertical
% sf = r1D; rf = -ones(size(r1D)); wf = w1D;

% % 1D LIFT
% Vf = bern_basis_1D(N,-1); Vf = Vf(:);
% V1D = bern_basis_1D(N,r1D);
% M1D = (V1D'*diag(w1D)*V1D);
% LIFT1D = M1D\Vf


% rp1D = linspace(-1,1,50); [rp sp] = meshgrid(rp1D); rp = rp(:); sp = sp(:);
% Vp = bern_quad(N,rp,sp);
% for i = 1:(N+1)^2
%     vv = Vp(:,i);
%     clf
%     %     color_line3(rp,sp,vv,vv,'.')
%     plot(r,s,'.');hold on;plot(r(i),s(i),'ro')
%     pause
% end
% return

%% randomly warp points + geofacs
r1 = [-1 1 -1 1]'; s1 = [-1 -1 1 1]';

% ep = .25; VX = r1 + ep*randn(size(r1)); VY = s1 + ep*randn(size(r1));
VX = r1 + [0 0 0 0]';
VY = s1 + [0 1 0 0]';
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

% plot(VX([1 2 4 3 1]),VY([1 2 4 3 1]),'o-');hold on
% plot(x,y,'s')
% plot(xf,yf,'s')
% quiver(xf,yf,nx,ny)
% axis equal; return

% compute operators

[V] = bern_quad(N,r,s);
Vf = bern_quad(N,rf,sf);
VS = sem_basis(N,r,s);
VSf = sem_basis(N,rf,sf);

% bb lift
M = V'*diag(w.*J)*V;
Mf = Vf'*diag(wf.*sJ)*Vf;
LIFT = M\Mf;
LIFT(abs(LIFT)<1e-6) = 0;
LIFT

% sem (diag) lift
MS = VS'*diag(w.*J)*VS;
MSf = VSf'*diag(wf.*sJ)*VSf;
LIFTS = MS\MSf;
LIFTS(abs(LIFTS)<1e-6) = 0;

imagesc(LIFT)
% keyboard
return

function V = sem_basis(N,r,s)

r1D = JacobiGL(0,0,N);
[ss rr] = meshgrid(r1D); rr = rr(:); ss = ss(:);
sk = 1;
for i = 0:N
    for j = 0:N
        V(:,sk) = JacobiP(rr,0,0,i).*JacobiP(ss,0,0,j);
        Vq(:,sk) = JacobiP(r,0,0,i).*JacobiP(s,0,0,j);
        sk = sk + 1;
    end
end
V = Vq/V;
