% orthogonal basis on mapped pyramidal elements. semi nodal

function [V Vr Vs Vt] = pyr_basis(N,r,s,t)

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

Np = (N+1)*(N+2)*(2*N+3)/6;
V = zeros(length(r),Np);
Vr = V; Vs = V; Vt = V;

ind = 1;
for k = 0:N
    [bq aq wab] = QCubature2D(k);
    VDM = QVandermonde2D(k,aq,bq); 
    [Vab Va Vb] = QVandermonde2D(k,a,b);

    Vab = Vab/VDM; DVa = Va/VDM; DVb = Vb/VDM;
        
    CNk = (N+2)/(2^(2*k+2)*(2*k+3));
    pNk = JacobiP(c,2*k+3,0,N-k);    
    pc = ((1-c)/2).^k.*pNk;
    dpNk = GradJacobiP(c,2*k+3,0,N-k);    
    dpc = ((1-c)/2).^k.*dpNk - (k/2^k)*((1 - c).^(k - 1)).*pNk;
    
    for ij = 1:(k+1)^2
        scale = 1/sqrt(CNk*wab(ij)); % normalization
        V(:,ind) = scale*Vab(:,ij).*pc;
        
        Vr(:,ind) = scale*(DVa(:,ij).*dadr).*pc;
        Vs(:,ind) = scale*(DVb(:,ij).*dbds).*pc;
        Vt(:,ind) = (DVa(:,ij).*dadt + DVb(:,ij).*dbdt).*pc + ...
            Vab(:,ij).*dpc*dcdt;
        Vt(:,ind) = scale*Vt(:,ind); 
        ind = ind+1;
    end
end

