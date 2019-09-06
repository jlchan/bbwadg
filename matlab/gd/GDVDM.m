function [V VX] = GDVDM(N,K,r)

nGhost = (N-1)/2;
Ka = (K/2 + nGhost); % half span + ghost point
VX = -Ka:Ka;

Npts = length(r);

% build ghost VDM
Np = length(-K/2:K/2);
V = zeros(Npts,Np);
for j = VX
    jid = Ka+1+j;
    for i = 1:Npts
        V(i,jid) = phi(r(i)-VX(jid),N);
    end
end

% modify extrapolation
ec = getExtrapCoeffs( N+1 );
for ig = nGhost:-1:1
    ileft  = 1+nGhost;
    iGhost = ileft-ig;
    for k = 2:length(ec)
        V(:,iGhost+k-1) = V(:,iGhost+k-1) + ec(k)*V(:,iGhost);
    end
end
for ig = nGhost:-1:1
    iright = K + nGhost + 1;
    iGhost = iright+ig;
    for k = 2:length(ec)
        V(:,iGhost-k+1) = V(:,iGhost-k+1) + ec(k)*V(:,iGhost);
    end
end
V = V(:,nGhost + (1:Np)); % extract non-ghost columns
VX = VX(nGhost+(1:Np));

end


function z = phi(xi,p)

z = 0;
a = -(p+1)/2;
b = -a;
if( xi >= a & xi < b )
    ia = p-ceil(xi-a)+1;
    ib = -(p-ceil(b-xi)+1);
    lp = [ia:-1:1,-1:-1:ib];
    
    den = 1;
    num = 1;
    for k = 1:length(lp)
        num = num*(xi+lp(k));
        den = den*lp(k);
    end
    z = num/den;
end

end