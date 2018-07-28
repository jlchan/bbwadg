function [V Vr Vs Vt] = Vandermonde3DHex(N,r,s,t)

V = zeros(length(r(:)),(N+1)^3);
sk = 1;
for i = 0:N
    for j = 0:N
        for k = 0:N
            V(:,sk) = JacobiP(r,0,0,i).*JacobiP(s,0,0,j).*JacobiP(t,0,0,k);
            Vr(:,sk) = GradJacobiP(r,0,0,i).*JacobiP(s,0,0,j).*JacobiP(t,0,0,k);
            Vs(:,sk) = JacobiP(r,0,0,i).*GradJacobiP(s,0,0,j).*JacobiP(t,0,0,k);
            Vt(:,sk) = JacobiP(r,0,0,i).*JacobiP(s,0,0,j).*GradJacobiP(t,0,0,k);
            sk = sk + 1;
        end
    end
end
