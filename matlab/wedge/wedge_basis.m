% orthogonal basis on mapped pyramidal elements.

function [V Vr Vs Vt] = wedge_basis(N,r,s,t)

r = r(:); s = s(:); t = t(:);

% convert to abc
a = 2*(r+1)./(1-t) - 1;
b = s;
c = t;

ids = abs(1-t)<1e-12;
a(ids) = -1; 

Np = (N+1)^2*(N+2)/2;
V = zeros(length(r),Np);
Vr = V; Vs = V; Vt = V;

r1D = JacobiGQ(0,0,N); invV1D = inv(Vandermonde1D(N,r1D));

ind = 1;
for j = 0:N
    
    pb = JacobiP(b,0,0,j);
    dpb = GradJacobiP(b,0,0,j);
    
    for i = 0:N
        scale = .5^(-(2*i+1)/2);        
        
        pa = JacobiP(a,0,0,i);
        dpa = GradJacobiP(a,0,0,i);
                        
        for k = 0:N-i
            
            pc = ((1-c)/2).^i.*JacobiP(c,2*i+1,0,k);
            V(:,ind) = scale*pa.*pb.*pc;
            
            % change of vars from rst to abc
            % dadr = 2/(1-t) = 2/(1-c)
            % dadt = 2*(r+1)/(1-t)^2 = (a+1)/(1-c) = (a+1)/2 * dadr
            % dp/dr = da/dr * dp/da
            % dp/dt = da/dt * dp/da
            
            % if i = 0, Pa = const, so dadr/s/t = 0
            dadr_pc = zeros(size(a)); dadt_pc = zeros(size(a));
            dpc = GradJacobiP(c,2*i+1,0,k);
            if i > 0
                dpc = (-i/2)*(((1-c)/2).^(i-1)).*JacobiP(c,2*i+1,0,k) + ...
                    ((1-c)/2).^i.*GradJacobiP(c,2*i+1,0,k);
                
                dadr_pc = ((1-c)/2).^(i-1).*JacobiP(c,2*i+1,0,k); % removable singularity
                dadt_pc = (a+1)/2.*dadr_pc;
            end
            
            Vr(:,ind) = scale*dpa.*pb.*dadr_pc;
            Vs(:,ind) = scale*pa.*dpb.*pc; % b = s
            Vt(:,ind) = scale*(dpa.*pb.*dadt_pc + pa.*pb.*dpc);
            
            ind = ind + 1;
        end
    end
end


