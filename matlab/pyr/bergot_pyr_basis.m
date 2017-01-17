function [V] = bergot_pyr_basis(N,r,s,t)

Np = (N+1)*(N+2)*(2*N+3)/6;

a = r./(1-t); b = s./(1-t);
ids = abs(t-1) < 1e-6;
a(ids) = 0; b(ids) = 0; % r,s = 0 for sym pyr at t = 1
c = 2*t-1;

% Vandermonde matrix
V = zeros(numel(r),Np);
off = 1;
for i = 0:N
    p1 = JacobiP(a,0,0,i);    
    for j = 0:N
        p2 = JacobiP(b,0,0,j);        
        mij = max(i,j);
        for k = 0:N-mij                        
            p3 = (1-t).^mij.*JacobiP(c,2*mij + 2,0,k);                                                
            V(:,off) = p1.*p2.*p3;
            
            off = off + 1;
        end
    end
end
