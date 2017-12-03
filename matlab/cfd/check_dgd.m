clear 

N = 6;
K = N+1;
dx = 1/(N+1);
rv = linspace(-1,1,N+2)'; % vertex
rc = linspace(-1+dx,1-dx,N+1)'; % center

x1 = [rv(1:end-1) rv(2:end)]';
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+10);
Nq = length(rq);
xq = Vandermonde1D(1,rq)/Vandermonde1D(1,[-1,1])*x1;

% hold on
% plot(rc,rc*0,'o')
% plot(rv,rv*0,'x')

N0 = N/2+1;
e = zeros(N+1,1); 
v0 = zeros(Nq,K);
for i = 1:N+1
    id = N+2-i;
    e(id)=1;
    v0(:,i) = Vandermonde1D(N,xq(:,N0))/Vandermonde1D(N,rc)*e;
    e(id)=0;
end
plot(xq,v0)
return

% bandwidth = N/2
Vqh = zeros(Nq*K,2*N+1);
for i = 1:2*N+1
    if i <= N
        vtmp = v0(:,end-i+1:end);        
        Vqh(1:(i)*Nq,i) = vtmp(:);
    elseif i==N+1
        Vqh(:,i) = v0(:);
    end
end
Vqh(:,N+2:end) = fliplr(flipud(Vqh(:,1:N)));

% V(:,1) = v0(:);
% for i = 2:N0
%     vtmp = v0(:,i:(N+1));
%     V(1:(N+1)*(K-i+1),i) = vtmp(:);
%     vtmp = v0(:,1:(N+2-i));
%     V(((i-1)*(N+1)+1):(N+1)*K,N0+i-1) = vtmp(:);
% end
    
fq = exp(-50*xq.^2);

h = 2/K;
wJq = repmat(wq(:),1,K); wJq = wJq(:)*h;
Pqh = (Vqh'*diag(wJq)*Vqh)\(Vqh'*diag(wJq));
plot(xq(:),Vqh*Pqh*fq(:),'x')

Vq = Vandermonde1D(N,rq)/Vandermonde1D(N,r);
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
plot(xq,Vq*Pq*fq,'o--')
hold on
plot(xq,fq,'-')
% plot(xq(:),sum(Vq,2),'o-')
%  for i = 1:size(Vq,2)
%     clf
%     plot(xq(:),Vq(:,i),'o')
%      pause
%  end

