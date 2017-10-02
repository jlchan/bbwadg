clear

N = 2;
[ra wa] = JacobiGQ(0,0,N);
[rb wb] = JacobiGQ(1,0,N);
[rc wc] = JacobiGQ(2,0,N);
[a b c] = meshgrid(ra, rb, rc);
[wa wb wc] = meshgrid(wa,wb,wc);
a = a(:); b = b(:); c = c(:);
wa = wa(:); wb = wb(:); wc = wc(:);
w = wa.*wb.*wc;
w = (4/3)*w/sum(w);

r = 0.5*(1+a).*.5.*(1-b).*(1-c)-1; 
s = 0.5*(1+b).*(1-c)-1; 
t = c;

% plot3(r,s,t,'o');return

V = bern_basis_tet(N,r,s,t);
M = V'*diag(w)*V;

VBhex = zeros(length(a),(N+1)^3);

[Eth Etw Etq] = get_Eth(N);

sk = 1;
for k = 0:N
    for j = 0:N
        for i = 0:N
            Va = bern_basis_1D(N,a);
            Vb = bern_basis_1D(N,b);
            Vc = bern_basis_1D(N,c);
            VBhex(:,sk) = Va(:,i+1).*Vb(:,j+1).*Vc(:,k+1);            
            sk = sk + 1;
        end
    end
end

norm(V-VBhex*Eth,'fro')
Mhex = VBhex'*diag(w)*VBhex; 

Pq = M\(Eth'*Mhex);
Pq(abs(Pq)<1e-8) = 0;
iEth = pinv(Eth);
iEth(abs(iEth)<1e-8) = 0;
U = Pq - iEth;
U(abs(U)<1e-8) = 0;

% sk = 1;
% rtri = 0.5*(1+a).*(1-b)-1;
% stri = b;
% Vtri = bern_basis_tri(N,rtri,stri);
% for k = 0:N
%     for i = 1:(N+1)*(N+2)/2;
%         VBwedge(:,sk) = Vtri(:,i).*Vc(:,k+1);
%         sk = sk + 1;
%     end
% end
% norm(V-VBwedge*Etw,'fro')
% Mwedge = VBwedge'*diag(w)*VBwedge;
% PqW = M\(Etw'*Mwedge);
% PqW(abs(PqW)<1e-8) = 0;
% iEtw = pinv(Etw);
% iEtw(abs(iEtw)<1e-8) = 0;
% UW = PqW - iEtw;
% UW(abs(UW)<1e-8) = 0;

a1D = JacobiGQ(0,0,N);
for i = 0:N
    E{i+1} = bern_basis_1D(N,a1D)\bern_basis_1D(N-i,a1D);
end

T = bern_basis_tet(N,r,s,t)\Vandermonde3D(N,r,s,t);
% T = bern_basis_tri(N,r,s)\Vandermonde2D(N,r,s);
% T1D = bern_basis_1D(N,a1D)\Vandermonde1D(N,a1D);
cond(Pq)
cond(T)
cond(M)

