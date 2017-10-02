clear
N = 4;

[rq sq tq w] = pyr_cubature(N);
[Vq Vr Vs Vt] = pyr_basis(N,rq,sq,tq);

% color_line3(a,b,c,tz,'.');%hold on; text(a+.1,b,c,num2str((1:length(a))'))

%%

f = @(r,s,t) r.*s.*t;
dfdr = @(r,s,t) s.*t;
dfds = @(r,s,t) r.*t;
dfdt = @(r,s,t) r.*s;

u = Vq'*(w.*f(rq,sq,tq)); % representational dofs for f
dudr = Vq'*(w.*dfdr(rq,sq,tq)); % representational dofs for f
duds = Vq'*(w.*dfds(rq,sq,tq)); 
dudt = Vq'*(w.*dfdt(rq,sq,tq)); 

ur = Vr'*(w.*f(rq,sq,tq));

%% define "semi-nodal" operators

[r s t] = pyr_semi_nodes(N);
[Dr Ds Dt] = pyr_deriv_ops(N);

norm(dudr-Dr*u)
norm(ur-Dr'*u)

norm(duds-Ds*u)
norm(dudt-Dt*u)

%% mapped

aa = .1;  ids = [1 3 4 2 5];
[VX VY VZ] = pyr_nodes(1); VX = VX(ids); VY = VY(ids); VZ = VZ(ids);
VX = VX + aa*randn(size(VX)); VY = VY + aa*randn(size(VX)); VZ = VZ + aa*randn(size(VX));
[xq,yq,zq,rxq,sxq,txq,~,~,~,~,~,~,Jq] = pyr_geom_factors(VX,VY,VZ,rq,sq,tq);

f = @(x,y,z) x.*y;
dfdx = @(x,y,z) y;

wJ = w.*Jq;
M = Vq'*diag(wJ)*Vq;
M(abs(M)<1e-8)=0;
% spy(M)
% keyboard

u = M\(Vq'*(wJ.*f(xq,yq,zq))); % representational dofs for f

dudx = Vq'*(wJ.*(rxq.*(Vr*u)+sxq.*(Vs*u)+txq.*(Vt*u))); %Vq'*(wJ.*dfdx(xq,yq,zq)); % representational dofs for f

uc = wJ.*(Vq*u);
Vxu = Vr'*(rxq.*uc) + Vs'*(sxq.*uc) + Vt'*(txq.*uc);

% Vxu = Vr'*(rx.*uc) + Vs'*(sx.*uc) + Vt'*(tx.*uc);

[x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = pyr_geom_factors(VX,VY,VZ,r,s,t);
% plot3(r,s,t,'.')

Sx = Vq'*diag(wJ)*(diag(rxq)*Vr + diag(sxq)*Vs + diag(txq)*Vt);
Dx = diag(rx)*Dr + diag(sx)*Ds + diag(tx)*Dt;
norm(M\Sx - Dx,'fro') % strong nodal operators

SVx = (Vr'*diag(rxq) + Vs'*diag(sxq) + Vt'*diag(txq))*diag(wJ)*Vq;
DVx = Dr'*diag(rx) + Ds'*diag(sx) + Dt'*diag(tx);

norm(M\SVx - diag(1./J)*DVx*diag(J),'fro') % weak nodal operators

