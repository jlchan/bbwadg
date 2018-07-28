clear -globals
Globals2D

a = .0625; % warping factor

N = 1;


% vol nodes
[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D); 
wq = wr(:).*ws(:);

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out

% face nodes
[rq1D wq1D] = JacobiGQ(0,0,N);
Nfq = length(rq1D);

e = ones(size(rq1D));
rfq = [rq1D; rq1D; -e; e]; 
sfq = [-e; e; rq1D; rq1D]; 

rfq = [rq1D(1) rq1D(1) rq1D(2) rq1D(2) -1 -1 1 1]';
sfq = [-1 1 -1 1 rq1D' rq1D']';
wfq = [wq1D; wq1D; wq1D; wq1D];
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vq1D = Vandermonde1D(N,rq1D)/V1D;

Vf1D = Vandermonde1D(N,[-1;1])/Vandermonde1D(N,rq1D);

Vf = Vandermonde2D(N,rfq,sfq)/Vandermonde2D(N,rq,sq);
