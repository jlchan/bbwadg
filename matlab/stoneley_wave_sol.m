syms k b1p b2p b1s b2s x y t B1 B2 B3 B4 w

u1a = (1i*k*B1*exp(-k*b1p*y) + k*b1s*B2*exp(-k*b1s*y))*exp(1i*(k*x-w*t));
u2a = (-k*b1p*B1*exp(-k*b1p*y) + 1i*k*B2*exp(-k*b1s*y))*exp(1i*(k*x-w*t));
u1b = (1i*k*B3*exp(k*b2p*y) - k*b2s*B4*exp(k*b2s*y))*exp(1i*(k*x-w*t));
u2b = (k*b2p*B3*exp(k*b2p*y) + 1i*k*B4*exp(k*b2s*y))*exp(1i*(k*x-w*t));

v1a = simplify(diff(u1a,t))
v2a = simplify(diff(u2a,t))
v1b = simplify(diff(u1b,t))
v2b = simplify(diff(u2b,t))

% return
% u1a = subs(u1a,t,0);
% u2a = subs(u2a,t,0);
% u1b = subs(u1b,t,0);
% u2b = subs(u2b,t,0);

u1ax = simplify(diff(u1a,x))
u2ay = simplify(diff(u2a,y))
u12axy = simplify(diff(u1a,y) + diff(u2a,x))

u1bx = simplify(diff(u1b,x))
u2by = simplify(diff(u2b,y))
u12bxy = simplify(diff(u1b,y) + diff(u2b,x))