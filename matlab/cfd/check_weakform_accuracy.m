num_discontin = 3;
a = [1 1 2 4];
c = [.5 .1 rand rand];
for N = 1:100
    [rq wq] = JacobiGQ(0,0,N);
    f = 0*rq;
    iex = 0;
    for i = 1:num_discontin
        f = f + a(i)*(rq > c(i));
        iex = iex + a(i)*(1-c(i));
    end
    val = wq'*f;
    h(N) = max(diff(rq));
    
    loglog(h(N),abs(val-iex),'o')
    hold on 
end
loglog(h,h,'--') % first order accuracy