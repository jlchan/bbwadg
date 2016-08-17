% clc

for N = 1:3
    [r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
%     [rq sq tq wq] = tet_cubature(2*N);
    [rq sq tq wq Fmaskq wfq] = tet_cubature_TW(N);
    
    Vq = Vandermonde3D(N,rq,sq,tq);
%     norm(Vq'*diag(wq)*Vq - eye(size(Vq,2)),'fro')
end
% plot3(rq,sq,tq,'o')
% plot3(rq(Fmaskq),sq(Fmaskq),tq(Fmaskq),'o')
% [r w] = JacobiGQ(0,0,N+2);
% [r w] = JacobiGR(1,0,N);
% JacobiP(r,1,0,N-1)'*(w.*JacobiP(r,1,0,N))


N = 2;
[rq sq tq wq Fmaskq wfq] = tet_cubature_TW(N);
rfq = rq(Fmaskq); sfq = sq(Fmaskq); tfq = tq(Fmaskq);
% [rfq sfq tfq wfq] = tet_surface_cubature(1);
plot3(rfq,sfq,tfq,'o')
Vfq = Vandermonde3D(N,rfq,sfq,tfq);
Msurf = Vfq'*diag(wfq)*Vfq;

