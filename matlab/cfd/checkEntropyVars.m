syms G kappa K
syms rho u E p

p = (G-1)*(E - rho*u^2/2);
eta = beta*(p / (rho^G))^(beta*(1-gamma));


%% numerical tests

clear

gamma = 1.4;
beta = 2;

rho = 2+rand; 
u = rand;
p = 1+rand;

m = rho*u;
E = rho*(p/(gamma-1) + u^2/2);

p = (gamma-1)*(E - rho*u^2/2);
rhoe = (E - rho*u^2/2);

pstar = -(p/(rho^gamma))^(1/(beta*(1-gamma)));
W = (pstar/p) * [E - 2*rhoe - p*(1+beta); -rho*u; rho];

Cgb = (2/(gamma-1) + (1+beta));
W(1) - pstar*(E/p - Cgb)
E - (W(1)/pstar+Cgb)*p

% check entropy vars = d(entropy)/dU
p = @(rho,m,E) (gamma-1)*(E - .5*m^2/rho);
h = @(rho,m,E) beta*(p(rho,m,E)/(rho^gamma))^(1/((1-gamma)*beta));
eta = @(rho,m,E) rho*h(rho,m,E);
delta = 1e-4;
W1 = (eta(rho+delta,m,E) - eta(rho-delta,m,E)) / (2*delta);
W2 = (eta(rho,m+delta,E) - eta(rho,m-delta,E)) / (2*delta);
W3 = (eta(rho,m,E+delta) - eta(rho,m,E-delta)) / (2*delta);

W./[W1;W2;W3]



%%
% clear
syms gamma
syms rho m E

u = m/rho;
p = (gamma-1)*(E - rho.*u.^2/2);

% flux Jacobian
F = [rho*u, rho*u^2 + p, (E+p)*u];
dFdU = [simplify(expand(diff(F,rho)));
    simplify(expand(diff(F,m)));
    simplify(expand(diff(F,E)))];

pstar = -(rho*p).^(-gamma/(gamma+1));
W = [pstar*E; -pstar*m; pstar*rho];

steps = 5;
dWdU = [transpose(simplify(expand(diff(W,rho)),steps));
     transpose(simplify(expand(diff(W,m)),steps));
     transpose(simplify(expand(diff(W,E)),steps))];

dUdW = simplify(expand(inv(dWdU)));
syms W1 W2 W3
rhop = ((gamma-1)*(W1.*W3 - W2.^2/2)).^((1+gamma)/(1-gamma));
pstar = -(rhop).^(gamma/(gamma+1));
U = [pstar.*W3; -pstar.*W2; pstar.*W1];

UU = matlabFunction(U);
% dUdW2 = [transpose(simplify(expand(diff(U,W1)),steps));
%     transpose(simplify(expand(diff(U,W2)),steps));
%     transpose(simplify(expand(diff(U,W3)),steps))];
% dUdW2 = subs(dUdW2,W1,W(1));
% dUdW2 = subs(dUdW2,W2,W(2));
% dUdW2 = subs(dUdW2,W3,W(3));

dFdW = simplify(dUdW*dFdU,steps);

e = simplify(dFdW-transpose(dFdW));
e
% e = subs(e,rho,1);
% e = subs(e,m,1);
% e = subs(e,E,1);
% e = subs(e,C,1/(beta*(1-gamma)));
% e = subs(e,gamma,1.4);
% e = subs(e,beta,2);

dFdW = subs(dFdW,rho,U(1));
dFdW = subs(dFdW,m,U(2));
dFdW = subs(dFdW,E,U(3));

dfdW11 = matlabFunction(dFdW(1,1))
dfdW12 = matlabFunction(dFdW(1,2))
dfdW13 = matlabFunction(dFdW(1,3))
dfdW22 = matlabFunction(dFdW(2,2))
dfdW23 = matlabFunction(dFdW(2,3))
dfdW33 = matlabFunction(dFdW(3,3))

% check homogeneous entropy property - dfdx = 1/(1+beta) * ((f_w *w)_x + 


%% test symmetrization with entropy vars

clear
syms rho m E 
syms gamma beta
assume(gamma,'real')
assume(gamma > 1)

u = m/rho;
p = (gamma-1)*(E - rho.*u.^2/2);
F = [rho*u, rho*u^2 + p, (E+p)*u];
dFdU = [simplify(expand(diff(F,rho)));
    simplify(expand(diff(F,m)));
    simplify(expand(diff(F,E)))];

% hughes entropy
H = -rho*log(p/rho^gamma);
% beta = 1;
H = -beta*rho*(p/rho^gamma)^(1/(beta*(1-gamma)));

W = [simplify(diff(H,rho));
    simplify(diff(H,m));
    simplify(diff(H,E))];

dWdU = [diff(W,rho) diff(W,m), diff(W,E)];
dUdW = inv(dWdU);

dFdW = dFdU*dWdU;
e = simplify(dFdW-transpose(dFdW));
% simplify(subs(e,gamma,1.4))
% dFdW = dFdU*dUdW;
matlabFunction 
matlabFunction(dWdU)
matlabFunction(dUdW)
matlabFunction(dFdW)


