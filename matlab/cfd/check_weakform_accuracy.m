clear
a = [1 1 2 4];
c = [.5 .1 rand];
c = [.5];
num_discontin = length(c);
for N = 1:250
    [rq wq] = JacobiGQ(0,0,N);
    f = 0*rq;
    iex = 0;
    for i = 1:num_discontin
        f = f + a(i)*(rq > c(i));
        iex = iex + a(i)*(1-c(i));
    end
    val = wq'*f;
    h(N) = max(diff(rq));
    
    err(N) = abs(val-iex);    
end
loglog(h,err,'o--') 
hold on 
loglog(h,h,'--') % first order accuracy

%%

N = 49;
[rq wq] = JacobiGQ(0,0,N);

rp = linspace(-1,1,1500)';
Vp = Vandermonde1D(N,rp)/Vandermonde1D(N,rq);

f = @(x) (x > 0);

figure(1)
% plot(rp,Vp*f(rq))
% hold on
plot(rq,f(rq).*wq,'ko')

% figure(2)
% plot(rq,wq,'o')
% hold on

%%
clear 
for N = 1:150
    [rq wq] = JacobiGQ(0,0,N);
    h(N) = max(diff(rq));
%     w(N) = mean(wq);
    err(1) = abs(rq(1)-wq(1)/2 + 1);
    for i = 1:N
        err(i+1) = abs((rq(i)+wq(i)/2) - (rq(i+1)-wq(i+1)/2));
    end
    err(N+2) = abs(rq(end)+wq(end)/2-1);    
    wmin(N) = min(err); 
    wsum(N) = sum(err);
end
loglog(h,wsum,'o')
hold on
loglog(h,wmin,'x')
loglog(h,h,'--')
loglog(h,h.^2,'--')
loglog(h,h.^4,'--')


%%

N = 10;
[rq wq] = JacobiGQ(0,0,N);

plot([-1 1],[0 0],'k.','markersize',16)
hold on
for i = 1:length(rq)
    plot(rq(i),0*max(wq)*i/(N+1),'o')
    
    plot(rq(i)+[-wq(i) wq(i)]/2,0*[1,1]*max(wq)*i/(N+1),'+-','markersize',10)    
end
   
%% 

Nq = 9;

[rq sq wq] = Cubature2D(Nq);

plot([-1 1 -1 -1],[-1 -1 1 -1],'o-')
hold on
t = linspace(0,2*pi,150)';
x = cos(t); y = sin(t);

% pi*r^2 = wq(i)
for i = 1:length(rq)
    plot(rq(i),sq(i),'.','markersize',16)
    ri = sqrt(wq(i)/pi);
    plot(ri*x+rq(i),ri*y+sq(i),'k-')
end