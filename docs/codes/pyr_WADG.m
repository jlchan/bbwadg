clear

% mapped pyr
aa = .25;  ids = [1 3 4 2 5];
[VX VY VZ] = pyr_nodes(1); VX = VX(ids); VY = VY(ids); VZ = VZ(ids);
VX = VX + aa*randn(size(VX)); VY = VY + aa*randn(size(VX)); VZ = VZ + aa*randn(size(VX));
    
for N = 1:9
    
    [rq sq tq w] = pyr_cubature(N);
    [Vq Vr Vs Vt] = pyr_basis(N,rq,sq,tq);
    
    [xq,yq,zq,rxq,sxq,txq,~,~,~,~,~,~,Jq] = pyr_geom_factors(VX,VY,VZ,rq,sq,tq);
    
    f = @(x,y,z) exp(x + y + z);
    
    wJ = w.*Jq;
    M = Vq'*diag(wJ)*Vq;
    M(abs(M)<1e-8)=0;
    
    u1 = M\(Vq'*(wJ.*f(xq,yq,zq))); % true L2 projection
    
    Mhat = Vq'*diag(w)*Vq;
    MW = Vq'*diag(w./Jq)*Vq;
    MWADG = Mhat*(MW\Mhat);
    u2 = MWADG\(Vq'*(wJ.*f(xq,yq,zq))); % true L2 projection
    
    err1(N) = sqrt(wJ'*(Vq*u1 - f(xq,yq,zq)).^2);
    err2(N) = sqrt(wJ'*(Vq*u2 - f(xq,yq,zq)).^2);
    
end
semilogy(err1,'o-','DisplayName','L2 projection')
hold on
semilogy(err2,'x-','DisplayName','WADG projection')