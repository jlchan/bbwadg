function test_lsc_hex

N = 3;
Nq = N+1;
[rq1D w1D] = JacobiGQ(0,0,Nq);
[rq sq tq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:); tq = tq(:);
[wr ws wt] = meshgrid(w1D); w = wr(:).*ws(:).*wt(:);

for ref = 1:4
    K1D = 2^(ref-1);
    [VY VX VZ] = meshgrid(linspace(-1,1,K1D+1));
    VX = VX(:); VY = VY(:); VZ = VZ(:);
    
    EToV = zeros(K1D^3,8);
    sk = 1;
    for k = 0:K1D-1
        for j = 0:K1D-1
            for i = 0:K1D-1
                quad_ids = [1 2 K1D+3 K1D+2];
                hex_ids = [quad_ids quad_ids+(K1D+1)^2];
                EToV(sk,:) = hex_ids + i + j*(K1D+1) + k*(K1D+1)^2;
                sk = sk + 1;
            end
        end
    end
    
    % % arnold type mesh
    h = (.5.^(ref-1));
    nlev = 2^(ref-1);
    t = linspace(-1,1,nlev+1);
    for lev = 2:nlev
        ids = find(abs(VZ-t(lev))<1e-8);
        
        [VYlev yids] = sort(VY(ids));
        VXlev = reshape(VX(ids(yids)),nlev+1,nlev+1);
        xids = zeros(nlev+1);
        for j = 1:nlev+1
            [~,xidj] = sort(VXlev(:,j));
            xids(:,j) = xidj + (j-1)*(nlev+1);
        end
        xids = xids(:);
        p = ids(yids(xids));
        p = reshape(p,nlev+1,nlev+1);
        toggle = [1 0 -1 0]; toggle = toggle(mod(lev-2,4)+1); % switch between add, don't move, subtract
        
        for line = 1:nlev+1
            perturb = repmat([-1 1],1,floor((nlev+1)/2));
            if mod(length(VZ(p)),2)==1
                perturb = [perturb -perturb(end)];
            end
            if mod(line,2)==1
                VZ(p(:,line)) = VZ(p(:,line)) + toggle*h*perturb(:);
            else
                VZ(p(:,line)) = VZ(p(:,line)) - toggle*h*perturb(:);
            end
        end
    end
    
    Vq = hex_basis(N,rq,sq,tq); %/hex_basis(N,r,s,t);
    Mref = Vq'*diag(w)*Vq;
    invMref = inv(Mref);
    a = 1; b = .75; c = .5;
    f = @(x,y,z) exp(a*x + b*y + c*z);
    
    K = size(EToV,1);
    proj_err = 0.0;
    lscproj_err = 0.0;
    wproj_err = 0.0;    
    for e = 1:K
        v = EToV(e,:);
        %v = v([1 2 3 4 8 5 1 2 6 7 3 7 8 5 6]); v = v(:);
        %plot3(VX(v),VY(v),VZ(v),'o-','linewidth',2);hold on
                
        [x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = hex_geom_factors(VX(v),VY(v),VZ(v),rq,sq,tq);
        fq = f(x,y,z);
        
        % L2 projection
        M = Vq'*diag(w.*J)*Vq;
        b = (Vq'*(w.*J.*fq));
        uproj = M\b;        
        errK = sum(w.*J.*(Vq*uproj - fq).^2);
        proj_err = proj_err + errK;
        
        % LSC-DG
        ulsc = invMref * (Vq'*(w.*sqrt(J).*fq));
        errK = sum(w.*J.*((Vq*ulsc)./sqrt(J) - fq).^2);
        lscproj_err = lscproj_err + errK;         
        
        % weighted projection
        weight = 1; %sqrt(J);
        MinvJ = (Vq'*diag(w./weight)*Vq);
%         M = Mref*(MinvJ\Mref); 
        invM = invMref * MinvJ * invMref;
        uwproj = invM*b;
        errK = sum(w.*J.*((Vq*uwproj)./J - fq).^2);
        wproj_err = wproj_err + errK;
    end
    err1(ref) = sqrt(proj_err);
    err2(ref) = sqrt(lscproj_err);
    err3(ref) = sqrt(wproj_err);
    % plot3(VX,VY,VZ,'o')
    % text(VX,VY,VZ,num2str((1:length(VX(:)))'))
end

h = .5.^(1:length(err1)); h = h(:);
loglog(h,err1,'o-','linewidth',2)
hold on
loglog(h,err2,'s-','linewidth',2)
loglog(h,err3,'x-','linewidth',2)
legend('L2 projection','LSC-DG','Weighted projection')
set(gca,'fontsize',14)

len = 1;
h = h(end-len:end);
errl2save = err1(end-len:end);
errprojsave = err3(end-len:end);
fit1 = [log(h(:)) ones(size(h(:)))]\log(errl2save(:));
fit2 = [log(h(:)) ones(size(h(:)))]\log(errprojsave(:));
title(sprintf('L2 proj rate = %f vs projected mass = %f\n',fit1(1),fit2(1)))



function [V Vr Vs Vt] = hex_basis(N,r,s,t)

Np = (N+1)^3;
V = zeros(length(r(:)),Np);
id = 1;
for i = 0:N

    p1 = JacobiP(r,0,0,i);
    dp1 = GradJacobiP(r,0,0,i);
    
    for j = 0:N
        
        p2 = JacobiP(s,0,0,j);
        dp2 = GradJacobiP(s,0,0,j);
        
        for k = 0:N

            p3 = JacobiP(t,0,0,k);
            dp3 = GradJacobiP(t,0,0,k);
            V(:,id) = p1.*p2.*p3;
            
            Vr(:,id) = dp1.*p2.*p3;
            Vs(:,id) = p1.*dp2.*p3;
            Vt(:,id) = p1.*p2.*dp3;
            
            id = id + 1;
        end
    end
end



function [x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = hex_geom_factors(VX,VY,VZ,r,s,t)

% get vertex nodal bases
r1 = [    -1     1    1     -1    -1     1    1     -1]';
s1 = [    -1    -1     1     1    -1    -1     1     1]'; 
t1 = [    -1    -1    -1    -1     1     1     1     1]';
% plot3(r1,s1,t1,'.')
% text(r1,s1,t1,num2str((1:8)'))

V1 = hex_basis(1,r1,s1,t1); 
invV = inv(V1);

[V1 Dr1 Ds1 Dt1] = hex_basis(1,r,s,t); % eval map @ cubature
Interp = V1*invV; Dr = Dr1*invV; Ds = Ds1*invV; Dt = Dt1*invV;

x = Interp*VX; y = Interp*VY; z = Interp*VZ;

xr = Dr*VX; yr = Dr*VY; zr = Dr*VZ;
xs = Ds*VX; ys = Ds*VY; zs = Ds*VZ;
xt = Dt*VX; yt = Dt*VY; zt = Dt*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

