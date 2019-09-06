% % 1D
% clear
% 
% K = 500;
% xv = linspace(-1,1,K+1)';
% x = .5*(xv(1:end-1)+xv(2:end));
% dx = 2/K;
% 
% Vr = [x.^(0:4) sin(pi*x) cos(pi*x) exp(x)];
% 
% sk = 1;
% for i = 1:size(Vr,2)
%     for j = i:size(Vr,2)
%         Vtest(:,sk) = Vr(:,i).*Vr(:,j);
%         sk = sk + 1;
%     end    
% end
% [wq id] = get_empirical_cubature(Vtest,ones(size(Vtest,1),1)*dx,1e-4,K/2);
% plot(x(id),wq,'o')

% 2D

clear
K = 50;
FinalTime = 2;
xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
x1D = x;
dx = 2/K;

e = ones(K-1,1);
Q = diag(e,1)-diag(e,-1);
% S(1,:) = 0; S(:,1) = 0;
% S(end,:) = 0; S(:,end) = 0;
% S(1,end) = -1; S(end,1) = 1;
% S(1,2) = 1; S(2,1) = -1;
% S(K-1,K) = 1; S(K,K-1) = -1;
Q(1,1) = -1;
Q(K,K) = 1;
Q = .5*Q;
M = dx*diag([.5; ones(K-2,1); .5]);
w1D = diag(M);

[x y] = meshgrid(x);
x = x(:); y = y(:);

N = 10;
sk = 1;
Vr = zeros(length(x(:)),(N+1)*(N+2)/2);
for i = 0:N
    for j = 0:N-i
        Vr(:,sk) = x.^i.*y.^j;
        sk=sk+1;
    end
end

sk = 1;
Vmass = 0*Vr; 
for i = 1:size(Vr,2)
    for j = 1:size(Vr,2)
        Vmass(:,sk) = Vr(:,i).*Vr(:,j);        
        sk = sk + 1;
    end    
end

tol = 1e-4;


% reduce target space
[Vmass,Smass,~] = rsvd(Vmass,10*N);
smass = diag(Smass);
smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
Vmass = Vmass(:,smass_energy > tol); % in 1D
disp('Done with Vmass svd')

[wq2 idv] = get_empirical_cubature(Vmass,dx^2*ones(size(x(:))),tol,inf);

ee = ones(size(x1D));
xf = [x1D; x1D(end)*ee; x1D; -x1D(end)*ee];
yf = [-x1D(end)*ee; x1D; x1D(end)*ee; x1D];
nx = [0*ee; ee; 0*ee; -ee];
ny = [-ee; 0*ee; ee; 0*ee];
wf0 = [w1D;w1D;w1D;w1D];


Mr = Vr'*diag(dx^2)*Vr;
Pr = Mr\(Vr'*diag(dx^2));

Qx = kron(Q,sparse(M));
Qy = kron(sparse(M),Q);
%norm(Pr'*Vr'*(Qx+Qx')*Vr*Pr - E'*(diag(nx.*wf0))*E,'fro')
%norm(Vr'*(Qx+Qx')*Vr - Vf'*(diag(nx.*wf0))*Vf,'fro')

e = ones(size(Qx,1),1);
ef = ones(size(wf0));

sk = 1;
Vf = zeros(length(xf(:)),(N+1)*(N+2)/2);
for i = 0:N
    for j = 0:N-i
        Vf(:,sk) = xf.^i.*yf.^j;
        sk=sk+1;
    end
end

E = Vf*Pr;

norm(Vr'*(Qx')*e - Vf'*diag(nx)*wf0,'fro')
norm(e'*(Qy)*Vr - ef'*(diag(ny.*wf0))*Vf,'fro')


sk = 1;
for i = 1:size(Vf,2)
    for j = 1:size(Vf,2)
        Vftest(:,sk) = Vf(:,i).*Vf(:,j);
        sk = sk + 1;
    end    
end

[Vmass,Smass,~] = svd(Vftest,0);
smass = diag(Smass);
smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
Vftest = Vmass(:,smass_energy > 1e-13); % in 1D
disp('Done with Vfmass svd')

% constrained lsq
C = [Vf'*diag(nx);
    Vf'*diag(ny)];
d = [Vr'*(Qx')*e; Vr'*(Qy')*e];

% linprog version
f = ones(size(wf0));

% |A*w - b| < delta
% A*w  < b + delta
% -A*w + b < delta
delta = tol;
disp('running linprog')
tic;
[wf id] = get_EQP_nodes(Vftest,wf0,delta,C,d);
toc

norm(C(:,id)*wf-d)

plot(xf(id),yf(id),'x')
hold on

% [wf id] = get_empirical_cubature_constrained(Vftest,wf0,C,d,tol,1e4)
return

%% separate xy constraints - doesn't seem to be any more accurate

% constrained lsq
Cx = [Vf'*diag(nx)];
Cy = [Vf'*diag(ny)];
dx = [Vr'*(Qx')*e];
dy = [Vr'*(Qy')*e];

% linprog version
f = ones(size(wf0));

% |A*w - b| < delta
% A*w  < b + delta
% -A*w + b < delta
disp('running linprog in x')
tic;
[wfx idx] = get_EQP_nodes(Vftest,wf0,tol,Cx,dx);
toc

disp('running linprog in y')
tic;
[wfy idy] = get_EQP_nodes(Vftest,wf0,tol,Cy,dy);
toc


plot(xf(idx),yf(idx),'bo')
hold on
plot(xf(idy),yf(idy),'rx')
return

%
%%
vv = wq2;
h = color_line3(x(idv),y(idv),vv,vv,'.');
set(h,'markersize',16)
hold on
vv = wf;
h = color_line3(xf(id),yf(id),vv,vv,'o');
set(h,'linewidth',2,'markersize',8)
for i = 1:length(xv)
    plot(xv(i)*ones(size(xv)),xv,'k-')
    plot(xv,xv(i)*ones(size(xv)),'k-')
end
axis equal
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
axis off

norm(Vftest(id,:)'*wf-Vftest'*wf0) 
norm(C(:,id)*wf-d)
sum(Vf(id,:)'*(wf.*nx(id))-Vr'*(Qx')*e)
sum(Vf(id,:)'*(wf.*ny(id))-Vr'*(Qy')*e)

