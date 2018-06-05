N = 4;
[r s t] = Nodes3D(N); 
[r s t] = xyztorst(r,s,t);
Np = length(r);
Nq = 2*N;
[rq sq tq wq] = tet_cubature(Nq);
V = Vandermonde3D(N,r,s,t);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
M =(Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));


[Vr Vs Vt] = GradVandermonde3D(N,r,s,t);
Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;

[rqtri sqtri wqtri] = Cubature2D(Nq);

e = ones(size(rqtri));
rfq = [rqtri;  -e; rqtri; -(1+rqtri+sqtri)];
sfq = [sqtri; rqtri;  -e; rqtri];
tfq = [-e;  sqtri; sqtri; sqtri];
wfq = [wqtri;wqtri;wqtri;wqtri];

Nfp = (N+1)*(N+2)/2; rf = r(1:Nfp); sf = s(1:Nfp);
Vfqf = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,rf,sf);
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
Lq = M\(Vfq'*diag(wfq));
VfPq = Vfq*Pq;

z = zeros(size(rqtri));
nrJ = [z; -e; z;  e];
nsJ = [z;  z; -e; e];
ntJ = [-e; z; z ; e];

Drq = Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq;
Dsq = Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq;
Dtq = Vq*Dt*Pq - .5*Vq*Lq*diag(ntJ)*Vfq*Pq;

WN = diag([wq;wfq]);
DNr = [Drq .5*Vq*Lq*diag(nrJ);
    -.5*diag(nrJ)*VfPq .5*diag(nrJ)];
DNs = [Dsq .5*Vq*Lq*diag(nsJ);
    -.5*diag(nsJ)*VfPq .5*diag(nsJ)];
DNt = [Dtq .5*Vq*Lq*diag(ntJ);
    -.5*diag(ntJ)*VfPq .5*diag(ntJ)];
QNr = WN*DNr;
QNs = WN*DNs;
QNt = WN*DNt;

subplot(1,2,1)
spy(DNr)
subplot(1,2,2)
spy(Vq*Lq)

length(rfq)+length(rq)
Np
% plot3(rfq,sfq,tfq,'o')
% hold on
% quiver3(rfq,sfq,tfq,nrJ,nsJ,ntJ)
% Vq*Dr*Pq - .5*Vq*
