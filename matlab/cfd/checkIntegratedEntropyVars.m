clear
syms rho m E 
syms gamma; % = 1.4;
assume(gamma,'real')
assume(gamma > 1)
% gamma = 1.4;

u = m/rho;
p = (gamma-1)*(E - rho.*u.^2/2);
F = [rho*u; rho*u^2 + p; (E+p)*u];
dfdU = [simplify(expand(diff(F,rho))), simplify(expand(diff(F,m))), simplify(expand(diff(F,E)))];

% hughes entropy 
H = -rho*log(p/rho^gamma);

W = [simplify(diff(H,rho));
    simplify(diff(H,m));
    simplify(diff(H,E))];

dWdU = simplify([diff(W,rho) diff(W,m), diff(W,E)]);
dUdW = simplify(inv(dWdU));
dfdW = simplify(dfdU*dUdW);
e = simplify(dfdW - transpose(dfdW))

return

params.Vars = [rho,m,E,gamma];
% dfdW = matlabFunction(dfdW,params)
% dWdU = matlabFunction(dWdU,params)
%dUdW = matlabFunction(dUdW,params)
dWdU11 = matlabFunction(dWdU(1,1),params)
dWdU12 = matlabFunction(dWdU(1,2),params)
dWdU13 = matlabFunction(dWdU(1,3),params)
dWdU22 = matlabFunction(dWdU(2,2),params)
dWdU23 = matlabFunction(dWdU(2,3),params)
dWdU33 = matlabFunction(dWdU(3,3),params)

dfdW11 = matlabFunction(dfdW(1,1),params)
dfdW12 = matlabFunction(dfdW(1,2),params)
dfdW13 = matlabFunction(dfdW(1,3),params)
dfdW22 = matlabFunction(dfdW(2,2),params)
dfdW23 = matlabFunction(dfdW(2,3),params)
dfdW33 = matlabFunction(dfdW(3,3),params)

return

%% check conversion

clear
% gamma = 1.4;
syms gamma

syms rho m E
u = m/rho;
p = (gamma-1)*(E - rho.*u.^2/2);

% hughes entropy 
H = -rho*log(p/rho^gamma);

W = [simplify(diff(H,rho));
    simplify(diff(H,m));
    simplify(diff(H,E))];
% params.Vars = [rho,m,E]; W = matlabFunction(W,params);

syms W1 W2 W3
alpha = ((gamma-1)/((-W3)^gamma))^(1/(gamma-1))*exp((-gamma + W1 - W2*W2/(2*W3)) / (gamma-1));
U = [-alpha*W3;
    alpha*W2;
    alpha*(1 - W2*W2/(2*W3))];
U = simplify(U);

% params.Vars = [W1,W2,W3]; U = matlabFunction(U,params);

rho = U(1);
u = U(2)/U(1);
E = U(3);
p = (gamma-1)*(E - rho.*u.^2/2);

F = [rho*u; rho*u^2 + p; (E+p)*u];
F = simplify(F);
dfdW = [simplify(diff(F,W1)), simplify(diff(F,W2)), simplify(diff(F,W3))];
simplify(dfdW-transpose(dfdW))
params.Vars = [W1,W2,W3,gamma]; 
dfdW11 = matlabFunction(dfdW(1,1),params)
dfdW12 = matlabFunction(dfdW(1,2),params)
dfdW13 = matlabFunction(dfdW(1,3),params)
dfdW22 = matlabFunction(dfdW(2,2),params)
dfdW23 = matlabFunction(dfdW(2,3),params)
dfdW33 = matlabFunction(dfdW(3,3),params)

%%
Nq = 4;
[rq wq] = JacobiGQ(0,0,Nq);
rq = (1+rq)/2;
wq = wq/2;

W1 = 2.1;
W2 = 2;
W3 = 5;

dFdW = 0;
for i = 1:length(rq)
    dFdW = dFdW + dfdW(W1*rq(i),W2*rq(i),W3*rq(i))*wq(i)*rq(i);
end
dFdW


