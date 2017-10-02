pyr_mesh;
N = 5;
e = 1;
v = EToV(e,:);
f1 = v([1 2 5]);

[r s] = Cubature2D(2*N);
[r s t] = map_triangle(VX(f1),VY(f1),VZ(f1),r,s);
plot3(r,s,t,'.')

f5 = v([1 2 3 4])
[r s] = meshgrid(JacobiGQ(0,0,N));
r = r(:); s = s(:);
[r s t] = map_quad(VX(f5),VY(f5),VZ(f5),r,s);
hold on
plot3(r,s,t,'.')

%%
[r s t] = pyr_nodes(1);
invV = inv(pyr_basis(1,r,s,t));

[a b c] = meshgrid(linspace(-1,.99,30));
[r s t] = pyr_abctorst(a,b,c);
Vp = pyr_basis(1,r,s,t);
color_line3(r,s,t,Vp*invV(:,5),'.')
%%
% check analytic formula
a = 2*(1+r)./(1-t)-1;
b = 2*(1+s)./(1-t)-1;
c = t;

L1 = (1-a).*(1-b).*(1-c)/8;
L2 = (1-a).*(1+b).*(1-c)/8;
L3 = (1+a).*(1-b).*(1-c)/8;
L4 = (1+a).*(1+b).*(1-c)/8;
L5 = (1+c)/2;
color_line3(r,s,t,L3,'.');colorbar
%%
syms r s t
a = 2*(1+r)./(1-t)-1;
b = 2*(1+s)./(1-t)-1;
c = t;

L1 = (1-a).*(1-b).*(1-c)/8;
L2 = (1-a).*(1+b).*(1-c)/8;
L3 = (1+a).*(1-b).*(1-c)/8;
L4 = (1+a).*(1+b).*(1-c)/8;
L5 = (1+c)/2;

simplify(diff(L4,r))
simplify(diff(L4,s))
simplify(diff(L4,t))
%%
pretty(simplify(L1))
pretty(simplify(L2))
pretty(simplify(L3))
pretty(simplify(L4))
pretty(simplify(L5))