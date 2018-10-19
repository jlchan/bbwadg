N = 8;

d = 2;

[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
[wr ws] = meshgrid(wq1D);
rq = rq(:); sq = sq(:); wq = wr(:).*ws(:);

[rq1D wq1D] = JacobiGL(0,0,N);
e = ones(size(rq1D));
rf = [rq1D; e; -rq1D; -e]; 
sf = [-e; rq1D; e; -rq1D]; 
wf = [wq1D; wq1D; wq1D; wq1D];

Vq = Vandermonde2D(N,rq,sq);
Vf = Vandermonde2D(N,rf,sf);

% plot(rq,sq,'o')
% hold on
% plot(rf,sf,'x')

lam = eig(Vf'*diag(wf)*Vf,Vq'*diag(wq)*Vq);

CT = d*N*(N+1)/2;
% CT = d*(N+1)*(N+2)/2;

fprintf('Computed CT = %f, full quad CT = %f, est CT = %f\n',max(abs(lam)),CT,CT*(2+1/N)^(d-1))

%% triangle

N = 2;
[rq sq wq] = Cubature2D(2*N);

Vq = Vandermonde2D(N,rq,sq);
M = Vq'*diag(wq)*Vq;

[rq1D wq1D] = JacobiGL(0,0,N);
rf = [rq1D; -rq1D; -ones(size(rq1D))];
sf = [-ones(size(rq1D)); rq1D; -rq1D];
wf = [wq1D;wq1D;wq1D];
Vf = Vandermonde2D(N,rf,sf);

Mf = Vf'*diag(wf)*Vf;

max(abs(eig(Mf,M)))

