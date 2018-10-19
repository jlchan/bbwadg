%% quads

useQuads = 1; mypath

N = 1;

[rq1D_vol wq1D_vol] = JacobiGQ(0,0,N);
[rq1D_face wq1D_face] = JacobiGL(0,0,N);

[rq sq] = meshgrid(rq1D_vol);
[wr ws] = meshgrid(wq1D_vol);
rq = rq(:); 
sq = sq(:);
wq = wr(:).*ws(:);

e = ones(size(rq1D_face));
rf = [rq1D_face; e; -rq1D_face; -e]; 
sf = [-e; rq1D_face; e; -rq1D_face]; 
wf = [wq1D_face; wq1D_face; wq1D_face; wq1D_face];

Vq = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Vf = Vandermonde2D(N,rf,sf);

M = Vq'*diag(wq)*Vq;
Mf = Vf'*diag(wf)*Vf;
K = Vr'*diag(wq)*Vr + Vs'*diag(wq)*Vs;

CT = max(abs(eig(Mf,M)))
CM = max(abs(eig(K,M)))

%% tris

clear
useQuads = 0; mypath

N = 5;

[rq sq wq] = Cubature2D(2*N);
[rq1D_face wq1D_face] = JacobiGQ(0,0,N);

rf = [rq1D_face; -rq1D_face; -ones(size(rq1D_face))];
sf = [-ones(size(rq1D_face)); rq1D_face; -rq1D_face];
wf = [wq1D_face; wq1D_face; wq1D_face];

Vq = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Vf = Vandermonde2D(N,rf,sf);

M = Vq'*diag(wq)*Vq;
Mf = Vf'*diag(wf)*Vf;
K = Vr'*diag(wq)*Vr + Vs'*diag(wq)*Vs;

CT = max(abs(eig(Mf,M)))

CM = max(abs(eig(K,M)))

CT_GLL = CT*(2+1/N)^(1/2)