r = JacobiGQ(0,0,3);

a = .5;
N = 7;
d = bern_basis_1D(N,a*r)./bern_basis_1D(N,r);
plot(log(d'),'o--')

%%

[r s] = Cubature2D(10);
[x y] = rstoxy(r,s);
a = .5;
[r2 s2] = xytors(a*x,a*y);
[a b] = rstoab(r,s);
[a2 b2] = rstoab(r2,s2);
plot(r,s,'o')
hold on
plot(r2,s2,'x')
return

N = 4;
d = bern_basis_tri(N,r,s)./bern_basis_tri(N,r2,s2);
d = d';
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
plot3(re,se,log(d(:,1)),'o')