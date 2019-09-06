clear

N = 3;
Ka = (N+1)/2;
VX = -Ka:Ka;
K = length(VX)-1;
[rr ww] = JacobiGQ(0,0,N);

h = VX(2)-VX(1);
rq = []; wq = [];
for e = 1:K
    rq = [rq; (1+rr)/2 + VX(e)];
    wq = [wq; h/2*ww];
end

% for i = 1:length(rq)
%     phivec(i,1) = phi(rq(i),N);
% end
% moment = 0;
% for k = N-1
%     moment = abs(wq'*((rq.^k).*phivec));
% end
% moment
% return

syms x
for k = 2:2:N-1
    moment = 0;
    phi2 = sym(0);
    a = ceil(N/2);
    for j = -(a-1):0 % half of the points
        ellj = sym(1);
        denom = sym(1);
        for i = -(a-1):a
            if i~=j
                ellj = ellj * (x-i);
                denom = denom * (j-i);
            end
        end
        ellj = ellj/denom;
        moment = moment + int(ellj*(x-j)^k,0,1);        
        
        phi2 = phi2 + int(ellj*ellj,0,1);
        
        if 0
            aa = 2;
            Ix = [-N,N];            
            %ezplot(ellj*(x-j),[Ix,-aa,aa])            
            ezplot(val,[Ix,-aa,aa])
            hold on
            plot([j ],[0],'o')
            val = moment;
            title(sprintf('ell_%d, moment = %g',j,val));
            pause
        end
        
        fprintf('N = %d, k = %d: moment = %f, phi2 = %f\n',N,k,moment,phi2)
    end
end

% clf
% ezplot(moment,[0,1])

return




return

ee = zeros(size(moment,2),1);
ee(1:2:end) = 1;
msum = moment*ee;
plot(rr,msum/max(abs(msum(:))))
hold on
plot(rr,moment*(1-ee)/max(abs(msum(:))))

% Vandermonde1D(10,rr)\msum

%%

re = linspace(-N,N,N+1)';
rp = linspace(-1,1,200)';
Vp = Vandermonde1D(N,rp)/Vandermonde1D(N,re);
rp2 = .5*(rp+1);
plot(rp,Vp(:,2).*rp2.^2)
hold on
plot(rp,Vp(:,1).*(rp2+1).^2)
plot(rp,Vp(:,3).*(1-rp2).^2)
plot(rp,Vp(:,4).*(-2+rp2).^2)


