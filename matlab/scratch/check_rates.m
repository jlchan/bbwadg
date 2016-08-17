clear all

% Driver script for solving the 3D IPDG wave equation
Globals3D;

% Order of polymomials used for approximation
N = 2;

sk = 1;
for ref = 2:4

%     % GMSH meshes
    run(sprintf('cubeTetra%d.m',ref))
%     [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ);

    % Initialize solver and construct grid and metric
    StartUp3D;
    
    % cubature test
    [rq, sq, tq, w] = tet_cubature(N);
    Vq = Vandermonde3D(N, rq, sq, tq)/V;
    xq = Vq*x; 
    yq = Vq*y; 
    zq = Vq*z; % phys cubature nodes        
    
    JK = J(1,:);
    wJ = w(:)*JK;
        
    % check that all PN is in space - should all be zero
    idx = 1;
    for i = 0:N
        for j = 0:N-i
            for k = 0:N-i-j
                uex = @(x,y,z) x.^i.*y.^j.*z.^k;
                u = (Vq'*diag(w)*Vq)\(Vq'*(diag(w)*uex(xq,yq,zq))); % projection       
                pterr = abs(Vq*u-uex(xq,yq,zq));
                testerr = max(pterr(:)); 
                if testerr>1e-12
                    keyboard
                end
                idx = idx + 1;
            end
        end
    end
    
    %uex = @(x,y,z) cos(pi*x).*cos(pi*y).*cos(pi*z);
    uex = @(x,y,z) sin(pi*x);
%     uex = @(x,y,z) x; % ok
%     uex = @(x,y,z) cos(0.01*pi*x);
    
    %     u = uex(x,y,z); % interpolant
    u = (V*V')\(Vq'*(diag(w)*uex(xq,yq,zq))); % projection
%     u = (Vq'*diag(w)*Vq)\(Vq'*(diag(w)*uex(xq,yq,zq))); % projection
        
    L2e = wJ.*(Vq*u-uex(xq,yq,zq)).^2;        
    
     h(sk) = max(1./Fscale(:));
%     h(sk) = max(JK*(4/3))^(1/3);
    
    err(sk) = sqrt(sum(L2e(:)));
    sk = sk + 1;
end
%h = .5.^(1:length(err));
loglog(h,err,'o-');hold on

F = [ones(length(h),1), log(h(:))];
C = (F'*F)\(F'*log(err(:)));
title(sprintf('Rate = %f',C(2)));
