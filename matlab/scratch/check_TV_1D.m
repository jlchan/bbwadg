function check_TV_1D


N = 5;
[rq wq] = JacobiGQ(0,0,N);
re = linspace(-1,1,N+1)';
rp = linspace(-1,1,500)';
r = JacobiGL(0,0,N);
V = bern_basis_1D(N,r);
Vq = bern_basis_1D(N,rq);
Vp = bern_basis_1D(N,rp);

f = @(r) r > 0;
u = V\f(r); %(Vq'*diag(wq)*Vq)\(Vq'*(wq.*f(rq)));

hold on
plot(rp,f(rp),'-','linewidth',2)
plot(rp,Vp*u,'--','linewidth',2)
plot(re,u,'o--')
legend('Exact function','Interpolant','Control points')


function TV = TV1D(N,u)

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end
