clear
for N = 5:9
    
    r = JacobiGL(0,0,N);
    [rq wq] = JacobiGQ(0,0,N+10);
    re = linspace(-1,1,N+1)';
    rp = linspace(-1,1,50)';
    VB = bern_basis_1D(N,r);
    Vp = bern_basis_1D(N,rp);
    VBe = bern_basis_1D(N,re);
    Vq = bern_basis_1D(N,rq);
    M = Vq'*diag(wq)*Vq;
    
    f1 = @(x) sin(x);
    f2 = @(x) 1 + (x > 0);
    f3 = @(x) (x > .1).*(x-.1);
    
    uB1 = VB\f1(r);
    uB2 = VB\f2(r);
    uB3 = VB\f3(r);
%     uB1 = M\(Vq'*(wq.*f1(rq)));
%     uB2 = M\(Vq'*(wq.*f2(rq)));
%     uB3 = M\(Vq'*(wq.*f3(rq)));
    
plot(rp,Vp*uB2,'-');
hold on;plot(re,uB2,'o');     return
    
    % find K with large TV
    TV1 = 0; TV2 = 0; TV3 = 0;
    for i = 1:N
        TV1 = TV1 + abs(uB1(i,:) - uB1(i+1,:));
        TV2 = TV2 + abs(uB2(i,:) - uB2(i+1,:));
        TV3 = TV3 + abs(uB3(i,:) - uB3(i+1,:));
    end
%     TV = TV./(N*max(abs(uB(:))));    
    TVN1(N) = TV1;
    TVN2(N) = TV2;
    TVN3(N) = TV3;
    
end
N = 1:N;
semilogy(N,TVN1,'o--')
hold on
semilogy(N,TVN2,'x--')
% semilogy(N,TVN3,'s--')

semilogy(N,exp(N*2/3),'--') % 2N/3 appears to be rate of TV growth for discontinuous functions