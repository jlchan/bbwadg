clear
N = 7;

r = JacobiGL(0,0,N);
re = linspace(-1,1,N+1)';
[rq wq] = JacobiGQ(0,0,2*N+2);
rp = linspace(-1,1,150)';

V = bern_basis_1D(N,r);
Vp = bern_basis_1D(N,rp);
Vq = bern_basis_1D(N,rq);
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

M = N;
rr = JacobiGL(0,0,M);
reM = linspace(-1,1,M+1)';
Ve = bern_basis_1D(N,reM);
E = bern_basis_1D(M,rr)\bern_basis_1D(N,rr);
reM = linspace(-1,1,M+1)';

a = .25;
f = @(x) (x-a).*(x > a);
% f = @(x) x > a;
f = @(x) exp(-25*(x-a).^2);
f = @(x) exp(-2.5*(x-a).^2);
% f = @(x) sin(3*(x-a));
u = Pq*f(rq);
% u = V\f(r);

plot(rp,f(rp))
hold on
plot(rp,Vp*u,'--');
plot(re,u,'o--')
% plot(reM,E*u,'o--')

%uu = Ve*u-E*u;
uu = E*u;
% uu = u;
w = ones(length(uu),1);
reu = linspace(-1,1,length(uu))';
% w = 1+(1+reu)/2 .* (1-reu)/2;

for i = 1:length(uu)-1
    du(i) = uu(i+1).*w(i+1)-uu(i).*w(i);
    TV(i) = abs(du(i));    
end
% TV = max(abs(du));

tol = 2/N;
for i = 2:length(uu)-1
    flag1 = abs(uu(i,:)-uu(i-1,:)) > tol;
    flag2 = abs(uu(i+1,:)-uu(i,:)) > tol;
    TV2(i) = abs(sign(uu(i,:) - uu(i-1,:)).*flag1 - sign(uu(i+1,:)- uu(i,:)).*flag2); % sign variations at each node
end
% du = abs(du);
% for i = 1:length(du)-1    
%     d2u(i) = abs(du(i+1)-du(i));
% end
stem(linspace(-1,1,length(TV)),TV,'x')
% stem(linspace(-1,1,length(TV2)),TV2,'x')
% stem(linspace(-1,1,length(d2u)),d2u,'o')
title(sprintf('Total var = %g, max error = %g\n',max(abs(du)),max(abs(f(rp)-Vp*u))))


%% directional check
addpath('./BBtests');

N = 7;
[rq sq wq] = Cubature2D(2*N);
[rp sp] = EquiNodes2D(150); [rp sp] = xytors(rp,sp);
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Vp = bern_basis_tri(N,rp,sp);
Vq = bern_basis_tri(N,rq,sq);
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

f = @(x,y) (x > 2*y-.5)*1.0;

u = Pq*f(rq,sq);

plot3(re,se,u,'o');
vv = f(rp,sp);
color_line3(rp,sp,vv,vv,'.')
vv = Vp*u;
color_line3(rp,sp,vv,vv,'.')

off = 0;
TVx = 0;
TVy = 0;
for i = 0:N
    Ni = N-i;
    Npi = Ni+1;
    idsx = (1:Npi) + off;
    TVx = TVx + TV1D(N-i,u(idsx,:));
    
    offj = 0;
    idsy = [];
    for j = 0:N-i
        idsy(j+1) = i + offj + 1;
        offj = offj + (N-j+1);
    end
    TVy = TVy + TV1D(N-i,u(idsy,:));
    off = off + Npi;
end

hold on
a = atan2(TVx,TVy);
plot([-1,cos(a)],[-1,sin(a)],'bo-')


