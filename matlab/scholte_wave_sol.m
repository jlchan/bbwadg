clear
syms k b1p b2p b1s b2s x y t B1 B2 B3 B4 w

u1a = (1i*k*B1*exp(-k*b1p*y)*exp(1i*(k*x-w*t)));
u2a = (-k*b1p*B1*exp(-k*b1p*y)*exp(1i*(k*x-w*t)));
u1e = (1i*k*B2*exp(k*b2p*y) - k*b2s*B3*exp(k*b2s*y)) * exp(1i*(k*x-w*t));
u2e = (k*b2p*B2*exp(k*b2p*y) + 1i*k*B3*exp(k*b2s*y)) * exp(1i*(k*x-w*t));

v1a = simplify(diff(u1a,t))
v2a = simplify(diff(u2a,t))
v1b = simplify(diff(u1e,t))
v2b = simplify(diff(u2e,t))

u1ax = simplify(diff(u1a,x))
u2ay = simplify(diff(u2a,y))
u12axy = simplify(diff(u1a,y) + diff(u2a,x))

u1bx = simplify(diff(u1e,x))
u2by = simplify(diff(u2e,y))
u12bxy = simplify(diff(u1e,y) + diff(u2e,x))