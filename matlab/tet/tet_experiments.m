N = 1;

%%
[rq sq tq w] = wedge_cubature(N);

Vq = wedge_basis(N,rq,sq,tq);

[r s t] = wedge_nodes(N);
V = wedge_basis(N,r,s,t);
M = inv(V*V');

Vq = Vq*inv(V);
norm(M-Vq'*spdiag(w)*Vq,'fro')
%%
[rq sq tq w] = tet_cubature(N);

% rstw = TN2012_12;
% rq = rstw(:,1);sq = rstw(:,2); tq = rstw(:,3); w = rstw(:,4);
% [rq sq tq] = xyztorst(rq,sq,tq); w = (4/3)*w/sum(w);

Vq = tet_basis(N,rq,sq,tq);

[r s t] = Nodes3D(N);[r s t] = xyztorst(r,s,t);
V = tet_basis(N,r,s,t);
M = inv(V*V');

Vq = Vq*inv(V);
norm(M-Vq'*spdiag(w)*Vq,'fro')