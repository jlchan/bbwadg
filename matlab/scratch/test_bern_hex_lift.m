
N = 4;
[r1D w1D] = JacobiGQ(0,0,N); 
% [r1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,r1D); w1D = sum(inv(V*V'),2);

[VB VrB] = bern_basis_1D(N,r1D); D = VB\VrB;

[r s t] = meshgrid(r1D); r = r(:); s = s(:); t = t(:); 
[wr ws wt] = meshgrid(w1D); w = wr(:).*ws(:).*wt(:);
[V Vr Vs Vt] = bern_hex(N,r,s,t);


[rf sf] = meshgrid(r1D); rf = rf(:); sf = sf(:); tf = -ones(size(rf));
[wrf wsf] = meshgrid(w1D); wf = wrf(:).*wsf(:);


% plot3(r,s,t,'o')
% hold on
% plot3(rf,sf,tf,'r*')

%% mapping

[r1 s1 t1] = meshgrid([-1 1]);  r1 = r1(:); s1 = s1(:); t1 = t1(:);

ep = 0; 
VX = r1 + ep*randn(size(r1)); VY = s1 + ep*randn(size(r1)); VZ = t1 + ep*randn(size(r1));
% VX = r1 + [0 0 0 0]'; VY = s1 + [0 1 0 0]';

% eval at cubature points
[V1 V1r V1s V1t] = bern_hex(1,r,s,t);
x = V1*VX; y = V1*VY; z = V1*VZ;
xr = V1r*VX; xs = V1s*VX; xt = V1t*VX; 
yr = V1r*VY; ys = V1s*VY; yt = V1t*VY;
zr = V1r*VZ; zs = V1s*VZ; zt = V1t*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);

rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

% face cubature pts for normals/ sJ
[V1f V1rf V1sf V1tf] = bern_hex(1,rf,sf,tf);
xf = V1f*VX; yf = V1f*VY; zf = V1f*VZ; 

xrf = V1rf*VX; xsf = V1sf*VX; xtf = V1tf*VX; 
yrf = V1rf*VY; ysf = V1sf*VY; ytf = V1tf*VY;
zrf = V1rf*VZ; zsf = V1sf*VZ; ztf = V1tf*VZ;
Jf = xrf.*(ysf.*ztf-zsf.*ytf) - yrf.*(xsf.*ztf-zsf.*xtf) + zrf.*(xsf.*ytf-ysf.*xtf);

rxf =  (ysf.*ztf - zsf.*ytf)./Jf; ryf = -(xsf.*ztf - zsf.*xtf)./Jf; rzf = (xsf.*ytf - ysf.*xtf)./Jf;
sxf = -(yrf.*ztf - zrf.*ytf)./Jf; syf =  (xrf.*ztf - zrf.*xtf)./Jf; szf = -(xrf.*ytf - yrf.*xtf)./Jf;
txf =  (yrf.*zsf - zrf.*ysf)./Jf; tyf = -(xrf.*zsf - zrf.*xsf)./Jf; tzf = (xrf.*ysf - yrf.*xsf)./Jf;

% face 1
nx = -txf; ny = -tyf; nz = -tzf;
sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*Jf;

% compute lift matrix
% V = bern_hex(N,r,s,t);
% Vf = bern_hex(N,rf,sf,tf); 

[rL sL tL] = meshgrid(JacobiGL(0,0,N)); % GLL nodes
rL = rL(:); sL = sL(:); tL = tL(:);

sk = 1;
for i = 0:N
    for j = 0:N
        for k = 0:N
            VN(:,sk) = JacobiP(rL,0,0,i).*JacobiP(sL,0,0,j).*JacobiP(tL,0,0,k);
            V(:,sk)  = JacobiP(r,0,0,i).*JacobiP(s,0,0,j).*JacobiP(t,0,0,k);
            Vf(:,sk) = JacobiP(rf,0,0,i).*JacobiP(sf,0,0,j).*JacobiP(tf,0,0,k);
            sk = sk + 1;
        end
    end
end
V = V*inv(VN); % nodal
M = V'*diag(w.*J)*V;
   
Mf = Vf'*diag(wf.*sJ)*Vf;

LIFT = M\Mf;
LIFT(abs(LIFT)<1e-6) = 0; 
% LIFTf = LIFT(:,1:(N+1)^2);

imagesc(LIFT)

% figure
% % plot(VX([1 2 4 3 1]),VY([1 2 4 3 1]),'o-');hold on
% plot3(x,y,z,'o');hold on
% plot3(xf,yf,zf,'r*')
% quiver3(xf,yf,zf,nx,ny,nz)
% % axis equal


