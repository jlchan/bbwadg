% look at 2D map b/w Legendre + BB

clear

N = 5;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
Np = length(r);

V = Vandermonde2D(N,r,s);
T = bern_basis_tri(N,r,s)\V;
d = max(abs(T),[],1);
d = T(1,:);
D2 = diag(d(:).^2);
T2 = T*diag(1./d);

%[rq sq wq] = Cubature2D(2*N+1);
[sq rq wq aq bq] = tri_tp_cubature(N);
% plot(aq,bq,'o');return
Vq = bern_basis_tri(N,rq,sq);
M = Vq'*diag(wq)*Vq;
%norm(inv(M)-T2*D2*T2','fro')
norm(inv(M)-T2*D2*T2','fro')

r1D = JacobiGL(0,0,N);
V1D = bern_basis_1D(N,r1D);
[rq1D wq1D] = JacobiGQ(0,0,N);
Vq1D = bern_basis_1D(N,rq1D);
M1D = Vq1D'*diag(wq1D)*Vq1D;

Etq = kron(Vq1D,Vq1D)\Vq;
Etq(abs(Etq)<1e-8) = 0;
% spy(Etq)

VDMq = Vandermonde2D(N,rq,sq);
% norm(eye(Np)-VDMq'*diag(wq)*VDMq,'fro')
Pquad = kron(M1D,M1D) \ kron(Vq1D',Vq1D') * diag(wq)*Vq; % projection onto 
Pquad(abs(Pquad)<1e-8) = 0;
Ptq = M\(Etq'*kron(M1D,M1D));
Ptq(abs(Ptq)<1e-8) = 0;
% spy(Ptq)

for i = 0:N
    %Ei{i} = bern_basis_1D(N-i+1,r1D)\bern_basis_1D(N-i,r1D);
    Ei{i+1} = bern_basis_1D(N,r1D)\bern_basis_1D(N-i,r1D);
end
[Eth Etw Etq] = get_Eth(N);

E = bern_basis_tri(N,r,s)\bern_basis_tri(N-1,r,s);
E(abs(E)<1e-8) = 0;

Tq = Etq*T2;
Tq(abs(Tq)<1e-8) = 0;
T1D = bern_basis_1D(N,r1D)\Vandermonde1D(N,r1D);
T1D = T1D*diag(1./T1D(1,:));
Eq = kron(T1D,T1D)\Tq; 
Eq(abs(Eq)<1e-8) = 0;
% spy(Tq)

offi = 0; 
for i = 0:N    
    idi = (1:N-i+1) + offi;
    offi = offi + N-i+1;    
    offj = 0;
    for j = 0:N
        idj = (1:N-j+1) + offj;        
        Mij{i+1,j+1} = M(idi,idj);
        Tij{i+1,j+1} = T(idi,idj);
        offj = offj + N-j+1;
    end   
end

C = zeros(N+1);
CT = {};
for i = 1:N+1
    for j = 1:N+1
        A = Ei{i}'*M1D*Ei{j};
        B = Mij{i,j};
        rat = A./B;
        if norm(rat-mean(rat(:)),'fro')<1e-8            
            C(i,j) = mean(rat(:));
        else
            keyboard
        end
        
        A = Ei{i}'*Tij{1,j};
        B = Tij{i,j};
        A(abs(B)<1e-8) = nan;
        B(abs(B)<1e-8) = nan;
        rat = A./B;
        rat(isnan(rat)) = 0;
        CT{i,j} = mean(rat,1);
    end
end

Ptq = M\(Etq'*kron(M1D,M1D));
A = T\(Ptq*Tq); A(abs(A)<1e-8) = 0;
imagesc(log(abs(A)))


% kron(Vq',Vq')*W*Vq

%% 3D legendre 

clear

N = 7;
[r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);

V = Vandermonde3D(N,r,s,t);
T = bern_basis_tet(N,r,s,t)\V;
d = max(abs(T),[],1);
d = T(1,:);
D2 = diag(d(:).^2);
T2 = T*diag(1./d);

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = bern_basis_tet(N,rq,sq,tq);
M = Vq'*diag(wq)*Vq;

r1D = JacobiGL(0,0,N);
[rt st] = Nodes2D(N); [rt st] = xytors(rt,st);
for i = 0:N
    %Ei{i} = bern_basis_1D(N-i+1,r1D)\bern_basis_1D(N-i,r1D);
    Ei{i+1} = bern_basis_tri(N,rt,st)\bern_basis_tri(N-i,rt,st);
end
[Eth Etw Etq] = get_Eth(N);


Tq = Eth*T2;
T1D = bern_basis_1D(N,r1D)\Vandermonde1D(N,r1D);
T1D = T1D*diag(1./T1D(1,:));
Eq = kron(T1D,kron(T1D,T1D))\Tq; 
Eq(abs(Eq)<1e-8) = 0;
spy(Eq)

[rq1D wq1D] = JacobiGQ(0,0,N);
Vq1D = bern_basis_1D(N,rq1D);
M1D = Vq1D'*diag(wq1D)*Vq1D;

Ptq = M\(Eth'*kron(M1D,kron(M1D,M1D)));
Th = Eth*T;
imagesc(log(abs(T\(Ptq*Th))))

