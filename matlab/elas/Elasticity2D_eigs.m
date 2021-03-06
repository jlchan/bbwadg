function ElasticityPulse2D

% clear all, clear
clear -global *

Globals2D

K1D = 4;
N = 4;
c_flag = 0;
FinalTime = 4;

filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
VX = (VX+1)/2;
VY = (VY+1)/2;
% VX = 2*VX;

StartUp2D;

% BuildPeriodicMaps2D(1,1);

% BCType = zeros(size(EToV));
% for e = 1:K
%     BCType(e,EToE(e,:)==e) = 6;
% end
%
% BuildBCMaps2D;
% nref = 1;
% for ref = 1:nref
%     Refine2D(ones(size(EToV)));
%     StartUp2D;
%     BuildBCMaps2D;
% end


[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%%
global mapBL vmapBL t0

t0 = .0;
% t0 = .075;
t0 = .05;

% find x = 0 faces
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
xfb = x(vmapP(fbids));
bfaces = sum(abs(xfb),1)<1e-8;

fbids = fbids(:,bfaces);
fbids = fbids(:);

mapBL = fbids;
vmapBL = vmapP(mapBL);

%%
global Nfld mu lambda Vq Pq tau useWADG
global C11 C12 C13 C22 C23 C33 rho


Nfld = 5; %(u1,u2,sxx,syy,sxy)

rho = ones(size(x))*1;
mu = ones(size(x))*1;
lambda = ones(size(x))*1;

mu0 = mu(1);
lambda0 = lambda(1);

% C = [2*mu0+lambda0       lambda0       0
%      lambda0       2*mu0+lambda0       0
%           0       0                mu0/2];

% mu = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% lambda = 1 + .5*sin(2*pi*x).*sin(2*pi*y);

% ids =  mean(y) > .375 & mean(y) < .625 & mean(x) > .5 & mean(x) < .75;
% mu(:,ids) = 10*mu(:,ids);
% lambda(:,ids) = 10*lambda(:,ids);


useWADG = 1;
if useWADG
    mu = Vq*mu;
    lambda = Vq*lambda;
    rho = Vq*rho;
    
    %     a = .5; b = 4;
    %     mu = 1 + a*sin(b*pi*xq).*sin(b*pi*yq);
    %     lambda = 1 + a*sin(b*pi*xq).*sin(b*pi*yq);
    %     color_line3(xq,yq,mu,mu,'.');return
    
    
    x0 = .5; y0 = .5;
    sc = .5;
    C11 = sc*(cos(2*pi*xq) + (yq > y0).*cos(2*pi*yq));
    C12 = sc*(sin(3*pi*xq) + (xq > x0).*sin(3*pi*yq));
    C13 = sc*(sin(3*pi*xq) + (xq > x0).*sin(3*pi*yq));
    C22 = sc*(cos(2*pi*xq) + (yq > y0).*cos(2*pi*yq));
    C23 = sc*(sin(3*pi*xq) + (xq > x0).*sin(3*pi*yq));
    C33 = sc*(cos(2*pi*xq) + (yq > y0).*cos(2*pi*yq));
    Cmin = 10;
    for i = 1:length(xq(:))
        C = [C11(i) C12(i) C13(i);
            C12(i) C22(i) C23(i);
            C13(i) C23(i) C33(i)];
        Cmin = min(Cmin,min(eig(C)));
    end
    aa = .1;
    C11 = C11 - Cmin + aa;
    C22 = C22 - Cmin + aa;
    C33 = C33 - Cmin + aa;
    
    % isotropic
    C11 = 2*mu+lambda;
    C12 = lambda;
    C13 = lambda;
    C22 = 2*mu+lambda;
    C23 = lambda;
    C33 = 2*mu+lambda;
    
    % for paper
    dmin = 1/10; dmax = 1;
    dmin = 1e-5; dmax = 1;
    for i = 1:length(xq(:))
        D = diag((dmax-dmin)*rand(3,1) + dmin);
        T = orth(randn(3));
        C = T*D*T';
        C11(i) = C(1,1);
        C12(i) = C(1,2);
        C13(i) = C(1,3);
        C22(i) = C(2,2);
        C23(i) = C(2,3);
        C33(i) = C(3,3);
    end
    
    iC11 = zeros(size(xq)); iC12 = zeros(size(xq));
    iC13 = zeros(size(xq)); iC22 = zeros(size(xq));
    iC23 = zeros(size(xq)); iC33 = zeros(size(xq));
    
    Cnorm = zeros(size(xq));
    
    Cmax = 0;
    for i = 1:length(xq(:))
        C = [C11(i) C12(i) C13(i);
            C12(i) C22(i) C23(i);
            C13(i) C23(i) C33(i)];
        Cnorm(i) = norm(C);
        Cmax = max(Cmax,norm(C));
        iC = inv(C);
        
        iC11(i) = iC(1,1);
        iC12(i) = iC(1,2);
        iC13(i) = iC(1,3);
        iC22(i) = iC(2,2);
        iC23(i) = iC(2,3);
        iC33(i) = iC(3,3);
    end
    
end

for tau0 = 0;
    
    for fld = 1:5
        cc = sqrt(max(rho(:).*Cnorm(:)));
        tau{fld} = tau0*ones(size(x));
        if fld > 2
            %tau{fld} = tau0./(2*mu+lambda);
%             tau{fld} = tau0/cc*ones(size(x));
            tau{fld} = tau0*ones(size(x));
            %                 tau{fld} = tau0./Cnorm;
            %         tau{fld} = tau0./Cnorm;
            %         tau{fld} = tau0./max(Cnorm(:))*ones(size(x));
        end
    end
    
    t0 = 0;
    u = zeros(Nfld*Np*K,1);
    rhs = zeros(Nfld*Np*K,1);
    A = zeros(Nfld*Np*K);
    ids = 1:Np*K;
    for i = 1:Nfld*Np*K
        u(i) = 1;
        for fld = 1:Nfld
            U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
        end
        rU = ElasRHS2D(U,0);
        u(i) = 0;
        for fld = 1:Nfld
            rhs(ids + (fld-1)*Np*K) = rU{fld};
        end
        A(:,i) = rhs(:);
        if (mod(i,100)==0)
            disp(sprintf('on col %d out of %d\n',i,Np*K*Nfld))
        end
    end
    lam = eig(A);
    
%     plot(lam,'bo','markersize',8,'linewidth',2)
    hold on
    plot(lam,'rs','markersize',8,'linewidth',2)
    %     hold on
    %     title(sprintf('tau = %f',tau0))
    %     axis equal
    %     drawnow
        set(gca,'fontsize',14)        

        grid on
    
    %         normA(ii,jj) = max(abs(lam));
    max(abs(lam))        
    %legend({'$d_{\min} = 10^{-1}$','$d_{\min} = 10^{-5}$'},'Interpreter','latex')
    
        keyboard
end

% surf(tauvvec,tausvec,normA')
% xlabel('\tau_v','fontsize',15)
% ylabel('\tau_s','fontsize',15)
% keyboard
return

function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda Vq Pq tau useWADG
global C11 C12 C13 C22 C23 C33 rho

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u;
    us = Ds*u;
    Ux{fld} = rx.*ur + sx.*us;
    Uy{fld} = ry.*ur + sy.*us;
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% velocity fluxes
nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};

opt=3;
if opt==1 % traction BCs
    %     global mapBL vmapBL t0
    %     if time < t0
    %         dU{1}(mapBL) = 2*(sin(pi*time/t0)*(time<t0) - U{1}(vmapBL));
    %         dU{2}(mapBL) = -2*U{2}(vmapBL);
    %     else
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    %     end
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
elseif opt==3 % zero velocity
    dU{1}(mapB) = -2*U{1}(vmapB);
    dU{2}(mapB) = -2*U{2}(vmapB);
    
end

% stress fluxes
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

% evaluate central fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms - reapply An
fp{1} = nx.*fc{3} + ny.*fc{5};
fp{2} = nx.*fc{5} + ny.*fc{4};
fp{3} = fc{1}.*nx;
fp{4} = fc{2}.*ny;
fp{5} = fc{2}.*nx + fc{1}.*ny;


flux = cell(5,1);
for fld = 1:Nfld
    flux{fld} = zeros(Nfp*Nfaces,K);
    flux{fld}(:) = fc{fld}(:) + tau{fld}(vmapM).*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu/2]
if useWADG
    for fld = 1:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = Pq*(rr{1}./rho);
    rhs{2} = Pq*(rr{2}./rho);
    %     rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    %     rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    %     rhs{5} = Pq*((mu) .* rr{5});
    rhs{3} = Pq*(C11.*rr{3} + C12.*rr{4} + C13.*rr{5});
    rhs{4} = Pq*(C12.*rr{3} + C22.*rr{4} + C23.*rr{5});
    rhs{5} = Pq*(C13.*rr{3} + C23.*rr{4} + C33.*rr{5});
else
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

return;

