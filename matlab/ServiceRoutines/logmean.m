% computes (aL-aR) / (log(aL)-log(aR))

function val = logmean(aL,aR)

xi = aL./aR;
f = (xi-1)./(xi+1);
u = f.^2;

F = log(xi)./2./f;
ids = abs(u) < 1e-4; % arbitrary
F(ids) = 1 + u(ids)/3 + u(ids).*u(ids)/5 + u(ids).*u(ids).*u(ids)/7;

val = (aL+aR)./(2*F);