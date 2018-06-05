useQuads = 1; mypath

N = 7;
[rq1D wq1D] = JacobiGQ(0,0,N);
V1D = Vandermonde1D(N,rq1D);
D1D = GradVandermonde1D(N,rq1D)/V1D;
Vf1D = Vandermonde1D(N,[-1;1])/V1D;
L1D = diag(1./wq1D)*(Vf1D');

[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D);
wq = wr(:).*ws(:);
V = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Dr = Vr/V;
Ds = Vs/V;
M = diag(wq);
e = ones(size(rq1D));
% rfq = [rq1D;e;-rq1D;-e];
% sfq = [-e;rq1D;e;-rq1D];

rfq = [rq1D; -rq1D; e;-e];
sfq = [-e;e;rq1D;-rq1D];

rtmp = [rq1D rq1D]';
stmp = [-e e]';
wtmp = [wq1D wq1D]';
rfq = [rtmp(:);stmp(:)];
sfq = [stmp(:);rtmp(:)];
wfq = [wtmp(:);wtmp(:)];
nrhat = [0*stmp(:);stmp(:)];
nshat = [stmp(:);0*stmp(:)];
Vf = Vandermonde2D(N,rfq,sfq)/V;

Lq = M\(Vf'*diag(wfq));

Lrq = Lq*diag(nrhat)*Vf;
norm(kron(L1D*diag([-1;1])*Vf1D,eye(N+1))-Lrq,'fro')
Lsq = Lq*diag(nshat)*Vf;
norm(kron(eye(N+1),L1D*diag([-1;1])*Vf1D)-Lsq,'fro')


plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on;

plot(rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
text(rq+.025,sq,num2str((1:(N+1)^2)'))
plot(rfq,sfq,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
text(rfq+.025,sfq,num2str((1:4*(N+1))'))
quiver(rfq,sfq,nrhat,nshat)
axis off
axis equal
