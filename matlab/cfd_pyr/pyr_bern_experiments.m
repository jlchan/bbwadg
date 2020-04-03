clear
N = 4;

[r s t] = pyr_nodes(N);
[rq sq tq wq] = pyr_cubature(2*N);
[rp sp tp] = pyr_equi_nodes(2*N);

[rt st wt] = Cubature2D(2*N);
et = ones(size(rt));
zt = 0*et;
[rq1D wq1D] = JacobiGQ(0,0,N);
[rquad squad] = meshgrid(rq1D);
rquad = rquad(:); 
squad = squad(:);
equad = ones(size(rquad));
zquad = 0*equad;
[wr ws] = meshgrid(wq1D); 
wquad = wr(:).*ws(:);

rf = [ rt; -rt;  rt; -et;  rquad];
sf = [-et;  st; -st;  rt;  squad];
tf = [ st;  rt;  st;  st; -equad];
wf = [ wt;  wt;  wt;  wt;  wquad];

nrJ = [ zt; et; zt; -et;  zquad];
nsJ = [-et; zt; et;  zt;  zquad];
ntJ = [ zt; et; et;  zt; -equad];

% [Vq Vrq Vsq Vtq] = bern_basis_pyr(N,rq,sq,tq);
[Vq Vrq Vsq Vtq] = pyr_basis(N,rq,sq,tq);
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
Vf = pyr_basis(N,rf,sf,tf);
Ef = Vf*Pq;

Qr = Vq'*diag(wq)*Vrq;
Br = diag(nrJ.*wf);

Qs = Vq'*diag(wq)*Vsq;
Bs = diag(nsJ.*wf);

Qt = Vq'*diag(wq)*Vtq;
Bt = diag(ntJ.*wf);

norm(Qr+Qr' - Vf'*Br*Vf,'fro')
norm(Qs+Qs' - Vf'*Bs*Vf,'fro')
norm(Qt+Qt' - Vf'*Bt*Vf,'fro')



