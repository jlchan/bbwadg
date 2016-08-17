clear
clear -global *

Globals3D;

% Order of polymomials used for approximation
N = 5;

% a = 1; f = @(x,y,z) sin(a*pi*x).*sin(a*pi*y).*sin(a*pi*z) + cos(x.*z);
f = @(x,y,z) exp(x + 2*y + .5*z + .1);% + log(x.*y.*z + 4);
%filenames = {'cubeTetra1.msh','cubeTetra2.msh','cubeTetra3.msh','cubeTetra4.msh'}; %,'cubeTetra5.msh'};
% filenames = {'cube1.msh','cube2.msh','cube3.msh','cube4.msh','cube5.msh'};
% filenames = {'cubeRef1.msh','cubeRef2.msh','cubeRef3.msh','cubeRef4.msh'};
filenames = {'sphere48.msh','sphere385.msh','sphere1780.msh','sphere10087.msh'};
% filenames = {'cubeUnstruc1.msh','cubeUnstruc2.msh','cubeUnstruc3.msh'};
% filenames = {'tet2_perturb.msh','tet3_perturb.msh','tet4_perturb.msh','tet5_perturb.msh'};
nfiles = length(filenames);
% nfiles = 3;
edgenum = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];

sk = 1;
for fileid = 1:nfiles
    filename = filenames{fileid};
    [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
    
    %     VX = VX + 2*randn(size(VX))/K;
    %     VY = VY + 2*randn(size(VX))/K;
    %     VZ = VZ + 2*randn(size(VX))/K;
    
    % Initialize solver and construct grid and metric
    StartUp3D;
    
    % if tet*_perturb.msh
    if (strcmp(filename(5:end),'_perturb.msh'))
        for e = 1:K
            if (min(J(:,e))<1e-8)
                p = [1 2 4 3];
                EToV(e,:) = EToV(e,p);
            end
        end
        StartUp3D
    end
    
    hMax = 0;
    for e = 1:K
        vx = VX(EToV(e,:));
        vy = VY(EToV(e,:));
        vz = VZ(EToV(e,:));
        h_edge = sqrt((vx(edgenum(:,1)) - vx(edgenum(:,2))).^2 + ...
            (vy(edgenum(:,1)) - vy(edgenum(:,2))).^2 + ...
            (vz(edgenum(:,1)) - vz(edgenum(:,2))).^2);
        hK(e) = max(h_edge);
        hK(e) = mean(h_edge);
    end
    hMax = max(hK); % option 1: max edge length over all elements - overpredicts convergence rate
    %     hMax = max(max(J)./max(sJ)); % option 2: ratio of J to sJ - underpredicts convergence rate
    %     hMax = max(J(:))^(1/3); % option 3: cube root of J - underpredicts convergence rate
    
    [rq sq tq wq] = tet_cubature(2*N+1);
    Vq = Vandermonde3D(N,rq,sq,tq)/V;
    xq = Vq*x;
    yq = Vq*y;
    zq = Vq*z;
    
    u = f(x,y,z); % nodal interpolant
    %     u = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*f(xq,yq,zq)); % L2 projection
    
    h(sk) = hMax;
    dofs(sk) = Np*K;
    e2 = 0;
    for e = 1:K
        e2 = e2+ J(1,e)*wq(:)'*(f(xq(:,e),yq(:,e),zq(:,e)) - Vq*u(:,e)).^2;
    end
    err(sk) = sqrt(e2);
    Kvec(sk) = K;
    sk = sk + 1;
end

%%
% color_line3(xq,yq,zq,log(abs(f(xq,yq,zq)-Vq*u)),'.')
% % color_line3(xq,yq,zq,Vq*u,'.');
% return
%%

if 0
    loglog(dofs,err,'x-');hold on
    r = N/3+1/3;
    rate = dofs.^(-r); rate = rate*err(1)/rate(1);
    loglog(dofs,rate,'--');hold on
    r = N/4+1/3;
    rate = dofs.^(-r); rate = rate*err(1)/rate(1);
    loglog(dofs,rate,'--');hold on
    C = polyfit(log(dofs(2:end)),log(err(2:end)),1);
    title(sprintf('Rate = %f, optimal = %f, non-optimal = %f\n',C(1), (N+1)/3, N/4 + 1/3))
else
    loglog(h,err,'o-');hold on
    C = polyfit(log(h(2:end)),log(err(2:end)),1);
    title(sprintf('Rate = %f, optimal = %f, non-optimal = %f\n',C(1), N+1,.75*N+1)) % N/4 + 1/3))
end

