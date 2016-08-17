clear
clear -global *

Globals2D;

N = 8;

filenames={'Grid/Maxwell2D/Maxwell1.neu','Grid/Maxwell2D/Maxwell05.neu',...
    'Grid/Maxwell2D/Maxwell025.neu','Grid/Maxwell2D/Maxwell0125.neu',...
    'Grid/Maxwell2D/Maxwell00625.neu'};
f = @(x,y) exp(x + y);

edgenum = [1 2; 2 3; 3 1];
sk = 1;
for id = 1:length(filenames)

    filename=filenames(id);
   
%     [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
    %[Nv, VX, VY, K, EToV] = unif_tri_mesh(id);
    StartUp2D;
    
    hMax = 0;
    for e = 1:K
        vx = VX(EToV(e,:));
        vy = VY(EToV(e,:));
        h_edge = sqrt((vx(edgenum(:,1)) - vx(edgenum(:,2))).^2 + ...
            (vy(edgenum(:,1)) - vy(edgenum(:,2))).^2);
        hK(e) = max(h_edge);
    end
    %     hMax = max(hK); % option 1: max edge length over all elements - overpredicts convergence rate
    %     hMax = max(max(J)./max(sJ)); % option 2: ratio of J to sJ - underpredicts convergence rate
    hMax = max(J(:))^(1/3); % option 3: cube root of J - underpredicts convergence rate
    
    [rq sq wq] = Cubature2D(2*N+2);
    Vq = Vandermonde2D(N,rq,sq)/V;
    xq = Vq*x; yq = Vq*y;
    
    %     u = f(x,y,z); % nodal interpolant
    u = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*f(xq,yq)); % L2 projection
    
    h(sk) = hMax;
    dofs(sk) = Np*K;
    e2 = 0;
    for e = 1:K
        e2 = e2+ J(1,e)*wq(:)'*(f(xq(:,e),yq(:,e)) - Vq*u(:,e)).^2;
    end
    err(sk) = sqrt(e2);
    Kvec(sk) = K;
    sk = sk + 1;
end

loglog(dofs,err,'x-');hold on
C = polyfit(log(dofs(2:end)),log(err(2:end)),1);
% loglog(h,err,'o-');hold on
% C = polyfit(log(h(2:end)),log(err(2:end)),1);
title(sprintf('Rate = %f, expected = %f\n',C(1), (N+1)/2))
