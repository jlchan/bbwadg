% test duffy expansion of BB
clear
% clc

N = 2;
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
    E{i+1} = bern_basis_1D(N,JacobiGL(0,0,N))\bern_basis_1D(N-i,JacobiGL(0,0,N));
end


VB1 = bern_basis_1D(N,a);
for j = 0:N
    for i = 0:N-j
        VB2 = bern_basis_1D(N,b)*E{j+1};
        Vq2(:,sk) = VB1(:,j+1).*VB2(:,i+1);
        sk = sk + 1;
    end
end

norm(Vq-Vq2,'fro') 
% return

% try TP
u = rand(Np,1);
uq = Vq*u;

off = 0;
for i = 0:N
    ids = (1:N-i+1) + off; 
    uTP(:,i+1) = E{i+1}*u(ids);
    off = off + N-i+1;
end


[a wa] = JacobiGQ(0,0,N);
VB = bern_basis_1D(N,a);
M1D = VB'*diag(wa)*VB;

% interp to qpts
uq2 = ((VB*uTP)*VB')'; 

norm(uq-uq2(:),'fro')

uq2 = (uq2.^2).*reshape(w,N+1,N+1);

% integrate
fTP = VB'*(VB'*uq2)';
% fTP = (inv(M1D')*uq2')*inv(M1D);
f1 = zeros(size(Vq,2),1);
off = 0;
for i = 0:N
    ids = (1:N-i+1) + off; 
    f1(ids) = E{i+1}'*fTP(:,i+1);
    off = off + N-i+1;
end

f2 = Vq'*(w.*(Vq*u).^2);
% M = Vq'*diag(w)*Vq;
% f2 = M\f2;

norm(f1-f2,'fro')

%% 3d test

clear

N = 2; 
Np = (N+1)*(N+2)*(N+3)/6;

Nq = N; % increase Nq to deal with quadratic weight
[ra wa] = JacobiGQ(0,0,Nq);
[rb wb] = JacobiGQ(0,0,Nq);
[rc wc] = JacobiGQ(2,0,Nq);

[a b c] = meshgrid(ra, rb, rc);
[wa wb wc] = meshgrid(wa,wb,wc);
a = a(:); b = b(:); c = c(:);
wa = wa(:); wb = wb(:); wc = wc(:);
w = wa.*wb.*wc.*(1-b);
w = (4/3)*w/sum(w);
r = 0.5*(1+a).*.5.*(1-b).*(1-c)-1; 
s = 0.5*(1+b).*(1-c)-1; 
t = c;

Vq = Vandermonde3D(N,r,s,t);
norm(Vq'*diag(w)*Vq - eye(size(Vq,2)),'fro')
% return
Vq = bern_basis_tet(N,r,s,t);

% for i = 1:size(Vq,2)
%     color_line3(r,s,t,Vq(:,i),'.')
%     pause
% end
for i = 0:N
    E{i+1} = bern_basis_1D(N,JacobiGL(0,0,N))\bern_basis_1D(N-i,JacobiGL(0,0,N));
end

Vc = bern_basis_1D(N,c);
sk = 1;
for k = 0:N
    Vb = bern_basis_1D(N,b)*E{k+1};
    for j = 0:N-k
        Va = bern_basis_1D(N,a)*E{j+k+1};
        for i = 0:N-j-k
            Vq2(:,sk) = Va(:,i+1).*Vb(:,j+1).*Vc(:,k+1);
            sk = sk + 1;
        end
    end
end
norm(Vq-Vq2,'fro')

% check application
u = (1:Np)'; %randn(Np,1);
uq = Vq*u;

[rt st] = Nodes2D(N); [rt st] = xytors(rt,st);
for i = 0:N
    Etri{i+1} = bern_basis_tri(N,rt,st)\bern_basis_tri(N-i,rt,st);
end

sk = 0;
for i = 0:N
    Nptri = (N-i+1)*(N-i+2)/2;
    ids = (1:Nptri)+ sk;
    uTri{i+1} = Etri{i+1}*u(ids);
    sk = sk + Nptri;    
end

uTP = zeros((N+1)^2,(N+1));
for k = 0:N
    utmp = zeros(N+1);
    off = 0;
    for i = 0:N
        ids = (1:N-i+1) + off;
        utmp(:,i+1) = E{i+1}*uTri{k+1}(ids);
        off = off + N-i+1;
    end
    uTP(:,k+1) = utmp(:);
end
uTP

[a wa] = JacobiGQ(0,0,Nq);
VB = bern_basis_1D(N,a);
[c wc] = JacobiGQ(2,0,Nq);
VBc = bern_basis_1D(N,c);
for k = 0:N
    uTPtmp = reshape(uTP(:,k+1),N+1,N+1);
    uqtmp = ((VB*uTPtmp)*VB')'; 
    uq2(:,k+1) = uqtmp(:);
end
uq2 = (VBc*uq2')';
norm(uq-uq2(:),'fro')

% uq2
c2 = reshape((1:(N+1)^3),(N+1)^2,N+1);
cuq = uq2.*c2;

% V^T
for k = 0:N
    uTPtmp = reshape(cuq(:,k+1),N+1,N+1);
    uqtmp = ((VB'*uTPtmp)*VB)'; 
    uint(:,k+1) = uqtmp(:);
end
uint = (VBc'*uint')';
uint
Eth = get_Eth(N);
b = Eth'*uint(:);

M = Vq'*diag(w)*Vq;
Pu = M\b
% b(1:end-1)

% Vq'*cuq(:)



% uq2 = (uq2(:).^2).*w;
% fTP = (VB'*uq2')*VB;
% % fTP = (inv(M1D')*uq2')*inv(M1D);
% f1 = zeros(size(Vq,2),1);
% off = 0;
% for i = 0:N
%     ids = (1:N-i+1) + off; 
%     f1(ids) = E{i+1}'*fTP(:,i+1);
%     off = off + N-i+1;
% end
% 
% f2 = Vq'*(w.*(Vq*u).^2);
% % M = Vq'*diag(w)*Vq;
% % f2 = M\f2;
% 
% norm(f1-f2,'fro')


%%

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
