for N = 1:150
    [rq wq] = JacobiGQ(0,0,N);
    f = 1.0*(rq > pi/4);
    val = wq'*f;
    h(N) = max(diff(rq));
    loglog(h(N),abs(val-(1-pi/4)),'o')
    hold on 
end
loglog(h,h,'--') % first order accuracy