clear
Globals2D
for N = 1:5;
    
Kvec = [8 16 32 64];
h = 2./Kvec;
sk = 1;
for K1D = Kvec    
    
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D);
    
    StartUp2D;
    
    % plotting nodes
    [rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
    Vp = Vandermonde2D(N,rp,sp)/V;
    xp = Vp*x; yp = Vp*y;
    % plot(xp,yp,'o'); return
        
    Nq = 2*N;
    [rq sq wq] = Cubature2D(Nq); % integrate u*v*c
    [rq2 sq2 wq2] = Cubature2D(Nq+2); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    Vq2 = Vandermonde2D(N,rq2,sq2)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); % J's cancel out
    xq = Vq*x; yq = Vq*y;
    Jq = Vq*J;
    
    rxJ = rx.*J; sxJ = sx.*J;
    ryJ = ry.*J; syJ = sy.*J;
    rxJ = Vq*rxJ; sxJ = Vq*sxJ;
    ryJ = Vq*ryJ; syJ = Vq*syJ;
    J = Vq*J;
    
    VqPq = Vq*Pq;
    
    gamma = 1.4;
    
    rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
    s = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
    V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - s(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
    V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
    V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
    V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));
    
    sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
    rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
    U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
    U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
    U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
    U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));
        
    rho = 2 + exp(.5*(xq+yq)).*sin(pi*xq).*sin(pi*yq);
    rhou = sin(pi*xq).*sin(pi*yq);
    rhov = sin(pi*xq).*sin(pi*yq);
    E = 2 + .5*(rhou.^2 + rhov.^2)./rho;
    
    rho = VqPq*rho;
    rhou = VqPq*rhou;
    rhov = VqPq*rhov;
    E = VqPq*E;
    
    q1 = Vq2*Pq*V1(rho,rhou,rhov,E);
    q2 = Vq2*Pq*V2(rho,rhou,rhov,E);
    q3 = Vq2*Pq*V3(rho,rhou,rhov,E);
    q4 = Vq2*Pq*V4(rho,rhou,rhov,E);
    
    rhoV = U1(q1,q2,q3,q4);
    rhouV = U2(q1,q2,q3,q4);
    rhovV = U3(q1,q2,q3,q4);
    EV = U4(q1,q2,q3,q4);
    
    rho = Vq2*Pq*rho;
    rhou = Vq2*Pq*rhou;
    rhov = Vq2*Pq*rhov;
    E = Vq2*Pq*E;
%     vv = rhovV;
%     color_line3(xq,yq,vv,vv,'.')
%     hold on
%     vv = rhov;
%     color_line3(xq,yq,vv,vv,'o')
    
    wJq = diag(wq2)*(Vq2*Pq*Jq);
    err2 = (rhoV - rho).^2 + (rhouV - rhou).^2  + (rhovV - rhov).^2  + (EV - E).^2;
    L2err(sk) = sqrt(sum(sum(wJq.*err2)));   
    sk = sk + 1;    
end
% disp(sprintf('N = %d',N))
print_pgf_coordinates(h,L2err)
C = [h(:).^0 log(h(:))]\log(L2err(:));
C(2)
C2 = diff(log(L2err(end-1:end)))/log(.5)
end
return
loglog(h,L2err,'bo--')
hold on
loglog(h,h.^(N+1),'rs-')
C = [h(:).^0 log(h(:))]\log(L2err(:));
C(2)
C2 = diff(log(L2err(end-1:end)))/log(.5)
% loglog(h,h.^(2*N+1),'k-')
