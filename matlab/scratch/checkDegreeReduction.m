% testing degree raising/lower
N = 3; i = 1;
[r s] = Nodes2D(N); [r s] = xytors(r,s);

V1 = bern_basis_tri(N,r,s);
V2 = bern_basis_tri(N-i,r,s);
E = V1\V2;

VDM1 = Vandermonde2D(N,r,s);
% VDM1 = VDM1*diag(1:size(V1,2));
T1 = inv(bern_basis_tri(N,r,s))*VDM1;
T1(abs(T1)<1e-8) = 0; 
[r s] = Nodes2D(N-i); [r s] = xytors(r,s);
T2 = inv(bern_basis_tri(N-i,r,s))*Vandermonde2D(N-i,r,s);

D = inv(T2)*E'*T1; D(abs(D)<1e-8) = 0;

[rq sq w] = Cubature2D(2*N);
V1q = bern_basis_tri(N,rq,sq); 
M1 = V1q'*diag(w)*V1q;

T = cell(N+1,1);
for ii = 0:N
    if ii<N
        r1D = JacobiGL(0,0,N-ii);
    else
        r1D = 0;
    end
    T{ii+1} = inv(bern_basis_1D(N-ii,r1D))*Vandermonde1D(N-ii,r1D);
end

%% interpolation for nodal

[r s] = Nodes2D(N); [r s] = xytors(r,s);
V1 = Vandermonde2D(N,r,s);
V2 = Vandermonde2D(N-i,r,s);
[r2 s2] = Nodes2D(N-i); [r2 s2] = xytors(r2,s2);

VDM2 = Vandermonde2D(N-i,r2,s2);
T = inv(VDM2)*inv(VDM2); % nodal to modal
T(abs(T)<1e-8) = 0;

V1*T*inv(VDM2)