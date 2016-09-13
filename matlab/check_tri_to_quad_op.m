% tet to hex op = sparse, approx (2*N+1) nnz per row up to N=9
% sum(c_i*Ei*Ei') * Etw'*Etq' * (quadrature on hex) * Eth
% - load Eth as sparse (can slab if needed)
% - do quadrature + project down
% - how to apply Etq', Etw' ?  

% look at tri to quad operator
clear
N = 6;
Np = (N+1)*(N+2)/2;
[r s w b a] = tri_tp_cubature(N);
Vq = Vandermonde2D(N,r,s);
% norm(Vq'*diag(w)*Vq-eye(size(Vq,2)),'fro')

% plot(r,s,'o'); text(a(:)+.1,b(:),num2str((1:length(a(:)))')); return

Vq = bern_basis_tri(N,r,s);

a = a(:); b = b(:);

sk = 1;
Vq2 = zeros(length(a),Np);

for i = 0:N
    E1D{i+1} = bern_basis_1D(N,JacobiGL(0,0,N))\bern_basis_1D(N-i,JacobiGL(0,0,N));
end

for i = 0:N
    Etmp = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
    Etmp(abs(Etmp)<1e-8) = 0;
    E{i+1} = Etmp;
end

% try TP - map tri to quad
u = zeros(Np,1);
for ii = 1:Np
    u(ii) = 1;
    
    uTP = zeros(N+1);
    off = 0;
    for i = 0:N
        ids = (1:N-i+1) + off;
        uTP(:,i+1) = E1D{i+1}*u(ids);
        off = off + N-i+1;
    end
    u(ii) = 0;
    Etq(:,ii) = uTP(:);
end
Etq(abs(Etq)<1e-8) = 0;

% spy(E{1}'*Etq')
% max(sum(abs(E{1}'*Etq')>0,2))

[a wa] = JacobiGQ(0,0,N);
VB = bern_basis_1D(N,a);
M1D = VB'*diag(wa)*VB;

% check  tet to wedge op
Np = (N+1)*(N+2)*(N+3)/6;

[r s] = Nodes2D(N);
[r s] = xytors(r,s);
for i = 0:N
    Etmp = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
    Etmp(abs(Etmp)<1e-8) = 0;
    E{i+1} = Etmp;
end

% tet to wedge
Etw = zeros((N+1)^2*(N+2)/2,Np);
u = zeros(Np,1);
for ii = 1:Np
    u(ii) = 1;
    sk = 0;
    for i = 0:N
        Nptri = (N-i+1)*(N-i+2)/2;
        ids = (1:Nptri)+ sk;
        uTri{i+1} = E{i+1}*u(ids);
        sk = sk + Nptri;
    end
    u(ii) = 0;
    uW = [uTri{:}];
    Etw(:,ii) = uW(:);
end

Eth = kron(eye(N+1),Etq)*Etw;



[r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
for i = 0:N
    Etmp = bern_basis_tet(N,r,s,t)\bern_basis_tet(N-i,r,s,t);
    Etmp(abs(Etmp)<1e-8) = 0;
    E{i+1} = Etmp;
end
j = 1;
hold on; 
spy(Eth,'.')

return


[rq sq tq wq ] = tet_cubature(2*N+1);
Vq = bern_basis_tet(N,rq,sq,tq);
M = Vq'*diag(wq)*Vq;
Mhex = kron(M1D,kron(M1D,M1D));
Pth = M\(Eth'*Mhex);

[r s w b a] = tri_tp_cubature(N);
Vq = bern_basis_tri(N,r,s);
Mtri = Vq'*diag(w)*Vq;
Mwedge = kron(M1D,Mtri);
Pth = M\(Etw'*Mwedge);
imagesc(Pth)

return

%% plot nnz per row
clear
for N = 1:12
    clear Etq Eth Etw
    Np = (N+1)*(N+2)/2;
    [r s w b a] = tri_tp_cubature(N);
    Vq = Vandermonde2D(N,r,s);
    % norm(Vq'*diag(w)*Vq-eye(size(Vq,2)),'fro')
    
    % plot(r,s,'o'); text(a(:)+.1,b(:),num2str((1:length(a(:)))')); return
    
    Vq = bern_basis_tri(N,r,s);
    
    a = a(:); b = b(:);
    
    sk = 1;
    Vq2 = zeros(length(a),Np);
    
    for i = 0:N
        E1D{i+1} = bern_basis_1D(N,JacobiGL(0,0,N))\bern_basis_1D(N-i,JacobiGL(0,0,N));
    end
    
    for i = 0:N
        Etmp = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
        Etmp(abs(Etmp)<1e-8) = 0;
        E{i+1} = Etmp;
    end
    
    % try TP - map tri to quad
    u = zeros(Np,1);
    for ii = 1:Np
        u(ii) = 1;
        
        uTP = zeros(N+1);
        off = 0;
        for i = 0:N
            ids = (1:N-i+1) + off;
            uTP(:,i+1) = E1D{i+1}*u(ids);
            off = off + N-i+1;
        end
        u(ii) = 0;
        Etq(:,ii) = uTP(:);
    end
    Etq(abs(Etq)<1e-8) = 0;
    
    % spy(E{1}'*Etq')
    % max(sum(abs(E{1}'*Etq')>0,2))
    
    [a wa] = JacobiGQ(0,0,N);
    VB = bern_basis_1D(N,a);
    M1D = VB'*diag(wa)*VB;
    
    
    % check  tet to wedge op
    % clear
    % N = 7;
    Np = (N+1)*(N+2)*(N+3)/6;
    
    [r s] = Nodes2D(N);
    [r s] = xytors(r,s);
    for i = 0:N
        Etmp = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
        Etmp(abs(Etmp)<1e-8) = 0;
        E{i+1} = Etmp;
    end
    
    % tet to wedge
    Etw = zeros((N+1)^2*(N+2)/2,Np);
    u = zeros(Np,1);
    for ii = 1:Np
        u(ii) = 1;
        sk = 0;
        for i = 0:N
            Nptri = (N-i+1)*(N-i+2)/2;
            ids = (1:Nptri)+ sk;
            uTri{i+1} = E{i+1}*u(ids);
            sk = sk + Nptri;
        end
        u(ii) = 0;
        uW = [uTri{:}];
        Etw(:,ii) = uW(:);
    end
    
    Eth = kron(eye(N+1),Etq)*Etw;
%     plot(N,max(sum(abs(Eth)>0,2)),'o')
    nnzEth(N) = max(sum(abs(Etw)>0,2));
%     hold on
end
N = 1:length(nnzEth);

Np = (N+1).*(N+2).*(N+3)/6;

plot(N,nnzEth,'x-')
hold on
plot(N,.15*N.^2+.75*N,'o--')
plot(N,2*N,'s--')
% plot(N,.2*N.^2,'--')

% spy(Eth)