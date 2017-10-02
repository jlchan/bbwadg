function [V Vr Vs Vt] = hex_basis(N,rin,sin,tin)

hybridgGlobalFlags % for useSEM flag

[r s t] = hex_nodes(N); % nvm, useSEM just used here
VDM = hex_VDM(N,r,s,t); 

% evaluate at rin,sin,tin
[V, Vr, Vs, Vt] = hex_VDM(N,rin,sin,tin); 
V = V/VDM; Vr = Vr/VDM; Vs = Vs/VDM; Vt = Vt/VDM;

function [V Vr Vs Vt] = hex_VDM(N,r,s,t)

Np = (N+1)^3;
V = zeros(length(r(:)),Np);
id = 1;
for i = 0:N

    p1 = JacobiP(r,0,0,i);
    dp1 = GradJacobiP(r,0,0,i);
    
    for j = 0:N
        
        p2 = JacobiP(s,0,0,j);
        dp2 = GradJacobiP(s,0,0,j);
        
        for k = 0:N

            p3 = JacobiP(t,0,0,k);
            dp3 = GradJacobiP(t,0,0,k);
            V(:,id) = p1.*p2.*p3;
            
            Vr(:,id) = dp1.*p2.*p3;
            Vs(:,id) = p1.*dp2.*p3;
            Vt(:,id) = p1.*p2.*dp3;
            
            id = id + 1;
        end
    end
end


