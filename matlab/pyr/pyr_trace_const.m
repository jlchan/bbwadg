clear
N = 5;

[rq sq tq wq] = pyr_cubature(N);
Np = (N+1)*(N+2)*(2*N+3)/6;

if 1
    % tri
    [rfq tfq wfq] = Cubature2D(2*N); sfq = -ones(size(rfq));
    Nfp = (N+1)*(N+2)/2;
    CN = 2*(N+1)*(N+2)/3 * sum(wfq)/sum(wq);
else
    % quad
    [rq1D wq1D] = JacobiGQ(0,0,N);
    [rfq sfq] = meshgrid(rq1D); rfq = rfq(:); sfq = sfq(:); tfq = -ones(size(rfq));
    [wr ws] = meshgrid(wq1D); wfq = wr(:).*ws(:);
    Nfp = (N+1)*(N+1);
    CN = (N+1)*(N+3)/3 * (sum(wfq)/sum(wq)); % quad;
end

[Vq] = pyr_basis(N,rq,sq,tq);
[Vfq] = pyr_basis(N,rfq,sfq,tfq);

% plot3(rfq,sfq,tfq,'o')
Mf = Vfq'*diag(wfq)*Vfq;

M = Vq'*diag(wq)*Vq;

max(eig(Mf,M))
(N+1) * Np/Nfp * (sum(wfq)/sum(wq))
CN




%%
N = 4;
[rq sq tq wq] = tet_cubature(2*N);
[rfq sfq wfq] = Cubature2D(2*N);
tfq = -ones(size(wfq));
Vq = Vandermonde3D(N,rq,sq,tq);
Vfq = Vandermonde3D(N,rfq,sfq,tfq);
M = Vq'*diag(wq)*Vq;
Mf = Vfq'*diag(wfq)*Vfq;

max(eig(Mf,M))
Nfp = (N+1)*(N+2)/2;
Np = (N+1)*(N+2)*(N+3)/6;
(N+1)* Np/Nfp * (sum(wfq)/sum(wq))