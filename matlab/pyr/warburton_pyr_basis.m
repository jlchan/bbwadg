function [V Vr Vs Vt] = warburton_pyr_basis(N,r,s,t)

Np = (N+1)*(N+2)*(2*N+3)/6;

% convert to abc
a = 2*(r+1)./(1-t)-1;
b = 2*(s+1)./(1-t)-1;
c = t;

ids = abs(1-t)<1e-12;
a(ids) = -1; b(ids) = -1;

% change of vars from a to b
dadr = 2./(1-t);
dbds = 2./(1-t);
dadt = (2*r + 2)./(1 - t).^2; 
dbdt = (2*s + 2)./(1 - t).^2; 
dcdt = 1;

% Vandermonde matrix
V = zeros(numel(r),Np);
Vr = zeros(numel(r),Np);
Vs = zeros(numel(r),Np);
Vt = zeros(numel(r),Np);

off = 1;
for i = 0:N
    p1 = JacobiP(a,0,0,i);    
    dp1 = GradJacobiP(a,0,0,i);
    for j = 0:N
        p2 = JacobiP(b,0,0,j);        
        dp2 = GradJacobiP(b,0,0,j);
        
        mij = (i+j);
        for k = 0:N-max(i,j) %mij                        
            p3 = (1-t).^mij.*JacobiP(c,2*mij + 2,0,k);
            
            dp31 = 0;                
            if mij > 0
                dp31 = -mij*(1-t).^(mij-1);
            end
            dp3 = dp31.*JacobiP(c,2*mij + 2,0,k) + (1-t).^mij.*GradJacobiP(c,2*mij + 2,0,k);
            
            V(:,off) = p1.*p2.*p3;
            
            va = dp1.*p2.*p3;
            vb = p1.*dp2.*p3;
            vc = p1.*p2.*dp3;
            
            Vr(:,off) = dadr.*va;
            Vs(:,off) = dbds.*vb;
            Vt(:,off) = dadt.*va + dbdt.*vb + dcdt.*vc;
            
            off = off + 1;
        end
    end
end
