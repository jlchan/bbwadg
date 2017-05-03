clear
N = 5;
[r w] = JacobiGQ(0,0,N);

alpha = 1; beta = 0;
for i = 0:N
    V(:,i+1) = JacobiP(r,alpha,0,i);
end

norm(V'*diag(w.*(1-r))*V - eye(N+1),'fro')

VB = bern_basis_1D(N,r);

T = VB\V;
x = 0:N;

we = @(alpha,beta,N) gamma((N+1) + alpha - (0:N)).*gamma(beta + 1+ (0:N))./(gamma(N+1- (0:N)).*gamma(1+ (0:N)));
W = diag(we(alpha,beta,N));
A = T'*W*T;
A(abs(A)<1e-8) = 0

% T'*A*T

%% tri connections

N = 3;

% [a b] = meshgrid(JacobiGQ(0,0,N));
% a = a(:); b = b(:);
% r = 0.5*(1+a).*(1-b)-1;
% s = b;
[r s w a b] = tri_tp_cubature(N);
% plot(r,s,'o')

a1D = JacobiGQ(0,0,N);

sk = 1;
V = zeros(length(a),(N+1)^2);
Va = zeros(N+1);
Vb = zeros(N+1);
for j = 0:N
    for i = 0:N
        V(:,sk) = JacobiP(a,0,0,i).*JacobiP(b,1,0,j);
        Va(:,i+1) = JacobiP(a1D,0,0,i);
        Vb(:,j+1) = JacobiP(a1D,1,0,j);
        sk = sk + 1;
    end
end
VB1D = bern_basis_1D(N,a1D);

VB = bern_quad(N,a,b) ;
Ta = VB1D\Va;
Tb = VB1D\Vb;

Vtri = Vandermonde2D(N,r,s);
VBtri = bern_basis_tri(N,r,s);
Etq = VB\VBtri; Etq(abs(Etq)<1e-8) = 0;  % bb tri to bb quad

T = VB\V; % convert from modal quad to bb quad

Tqtri = VB\Vtri; % convert from modal tri to bb quad
Ttri = VBtri\Vtri; % convert from modal tri to bb tri

norm(Etq*Ttri-Tqtri,'fro')
norm(T-kron(Tb,Ta),'fro')

we = @(alpha,beta,N) gamma((N+1) + alpha - (0:N)).*gamma(beta+1+ (0:N))./(gamma(N+1- (0:N)).*gamma(1+ (0:N)));
w1D = we(alpha,beta,N);
W = ones(N+1,1)*w1D(:)';
A = T'*diag(W(:))*T;
A(abs(A)<1e-8) = 0;
norm(A-diag(diag(A)),'fro')

%% try to apply inv(T) on quad

% [ae be] = meshgrid(linspace(-1,1,25));
% ae = ae(:); be = be(:);
% VB = bern_quad(N,ae,be);
% plot3(ae,be,VB(:,6),'.')
% return

clear

% M \ (BNi,f)

N = 2;
[r s w a b] = tri_tp_cubature(N);
% plot(a,b,'o')
% text(a,b,num2str((1:length(a))'))
% return

V = Vandermonde2D(N,r,s);
VB = bern_basis_tri(N,r,s);
M = VB'*diag(w)*VB;

a1D = JacobiGQ(0,0,N);

sk = 1;
Vquad = zeros(length(a),(N+1)^2);
Va = zeros(N+1); Vb = zeros(N+1);
for j = 0:N
    for i = 0:N
        Va(:,j+1) = JacobiP(a1D,0,0,j);
        Vb(:,j+1) = JacobiP(a1D,1,0,j);
        Vquad(:,sk) = JacobiP(a,0,0,i).*JacobiP(b,1,0,j);
        sk = sk + 1;
    end
end
VB1D = bern_basis_1D(N,a1D);
VBquad = bern_quad(N,a,b);
Ta = VB1D\Va;
Tb = VB1D\Vb;
Tquad = VBquad\Vquad; % bb to modal

alpha = 1; beta = 0;
we = @(alpha,beta,N) gamma((N+1) + alpha - (0:N)).*gamma(beta+1+ (0:N))./(gamma(N+1- (0:N)).*gamma(1+ (0:N)));
w1D = we(alpha,beta,N);
W = ones(N+1,1)*w1D(:)';

norm(Tquad-kron(Tb,Ta),'fro')

VB1D = bern_basis_1D(N,a1D);
VBquad = bern_quad(N,a,b);

Etq = VBquad\VB; Etq(abs(Etq)<1e-8) = 0;  % bb tri to bb quad

T = VB\V; % modal to bb
iT = V\VB; % bb to modal
D = T'*T;

Tq = Etq*T;

Np = size(V,2);
f = (1:Np)';
u = M\f;

% pinv(Etq)*Etq*T*pinv(Etq) * Etq*b
norm(inv(M)- T*T','fro')
norm(M-iT'*iT,'fro')

Pq = (M\(Etq'*VBquad'));
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
[be ae] = meshgrid(linspace(-1,1,N+1));
ae = ae(:); be = be(:);

return

Tq = Etq*T;
for i = 1:size(T,2)
    clf
    plot3(ae,be,Tq(:,i),'s--')
%     hold on
%     plot3(re,se,T(:,i),'o--')
    title(sprintf('id = %d\n',i))
    pause
end


