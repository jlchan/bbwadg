function Vp = Vandermonde2DQuad(N,xp,yp)

Vp = zeros(length(xp(:)), (N+1)*(N+1));
sk = 1;
for i=0:N    
    for j=0:N        
        Vp(:,sk) = JacobiP(xp(:), 0, 0, i).*JacobiP(yp(:), 0, 0, j);        
        sk = sk+1;        
    end    
end