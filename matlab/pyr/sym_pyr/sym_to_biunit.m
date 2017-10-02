% converts symmetric pyramid nodes to biunti nodes
function [r s t] = sym_to_biunit(r,s,t)

% biunit verts
rb = [-1 1 1 -1 -1]; sb = [-1 -1 1 1 -1]; tb = [-1 -1 -1 -1 1];
rb = rb(:); sb= sb(:); tb = tb(:);
brst = [rb sb tb];

% sym verts
rs = [-1 1 1 -1 0]; ss = [-1 -1 1 1 0]; ts = [0 0 0 0 1];
rs = rs(:); ss = ss(:); ts = ts(:);

% make vertex map
V1 = sym_pyr_basis(1,rs,ss,ts); % eval map @ cubature
map = V1\brst;

% new coordinates
V = sym_pyr_basis(1,r,s,t); % eval map @ cubature
r = V*map(:,1); s = V*map(:,2); t = V*map(:,3);

% bergot basis on the symmetric pyramid
function V = sym_pyr_basis(N,r,s,t)

Nfp = 3*N^2+2;
Nip = 1/6*(N-1)*(N-2)*(2*N-3);
Np = Nfp + Nip;

a = r./(1-t);
b = s./(1-t);
ids = abs(t-1)<1e-12;
a(ids) = 0; b(ids) = 0;
% Vandermonde matrix
V = zeros(numel(r),Np);
off = 1;
for i = 0:N
    p1 = JacobiP(a,0,0,i);    
    for j = 0:N
        p2 = JacobiP(b,0,0,j);        
        mij = max(i,j);
        for k = 0:N-mij                        
            if mij==0
                p3 = JacobiP(2*t-1,2*mij + 2,0,k);
            else
                p3 = (1-t).^mij.*JacobiP(2*t-1,2*mij + 2,0,k);
            end
            V(:,off) = p1.*p2.*p3;
            
            off = off + 1;
        end
    end
end
