% test_bern_tri_multiply

N = 3;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
V = bern_basis_tri(N,r,s);
V1 = bern_basis_tri(1,r,s);