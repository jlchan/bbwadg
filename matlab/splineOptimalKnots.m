clear
NN = 5;
KK = 32;

VX = linspace(-1,1,KK+1);
t = [VX(1)*ones(1,NN) VX VX(end)*ones(1,NN)]; % open knot vec

% interpolation (greville) points
if KK==1
    r = JacobiGL(0,0,NN);
else
    for i = 1:NN+KK
        r(i,1) = mean(t((i+1):(i+NN)));
    end
end
[rq wq] = JacobiGQ(0,0,NN);
rp = linspace(-1,1,25)';
rp = rp*.999;
% plot(r,r*0,'o')

h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,KK)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*KK,1);

rB = map(r); rB = rB(:);
rBq = map(rq); rBq = rBq(:);
wBq = repmat(wq,1,KK).*h(rq)/2; wBq = wBq(:);

R = bspline_basismatrix(NN+1,t,rB); % for local extraction to Lagrange dofs
R(abs(R)<1e-8) = 0;

% local basis
V = Vandermonde1D(NN,r);
Vq = Vandermonde1D(NN,rq)/V;
Dr = GradVandermonde1D(NN,r)/V;
Vp = Vandermonde1D(NN,rp)/V;

%% approx operators

NB = 4; % degree of spline space
Ksub = 16; % number of splines
N = NB+Ksub-1;

Bq = kron(speye(KK),Vq)*R;
%DBr = kron(spdiag(2./diff(VX)),Dr)^NB;
DBr = kron(speye(KK),Dr^NB);
Brq = kron(speye(KK),Vq)*DBr*R;

M = Bq'*diag(wBq)*Bq;
K = Brq'*diag(wBq)*Brq;

% rescale for numerical - assume equispaced knots for fine grid splines
K = K*(KK^(NB));
M = M/(KK^(NB));

% apply BCs
for i = 1:NB
    id = i;
    M(id,:) = 0; M(:,id) = 0;  M(id,id) = 1;
    K(id,:) = 0; K(:,id) = 0;  K(id,id) = 1;
    id = size(M,2)-i+1;
    M(id,:) = 0; M(:,id) = 0;  M(id,id) = 1;
    K(id,:) = 0; K(:,id) = 0;  K(id,id) = 1;
end
[W D] = eig(M,K);
[lam p] = sort(diag(D),'descend');
W = W(:,p);

VX = linspace(-1,1,Ksub+1);
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
re = linspace(-1,1,N+1)';
reKsub = linspace(-1,1,Ksub+1)';
for ii = 1:25
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)];
    Ve = bspline_basismatrix(NB+1, t, reKsub);
    VX = Ve*re;
    VX = VX(:)';
end

id = 2*NB+Ksub; % skip over BCs

% compute optimal knots as roots
VXK = linspace(-1,1,KK+1);
tK = [VXK(1)*ones(1,NN) VXK VXK(end)*ones(1,NN)]; % open knot vec
VXopt = zeros(size(VX));
for ii = 1:length(VX)
    VXopt(ii) = fzero(@(r) bspline_basismatrix(NN+1, tK, r)*W(:,id), VX(ii));
end


%% plotting and errors with smoothed knots

% compute errors
VX = linspace(-1,1,Ksub+1);
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
re = linspace(-1,1,N+1)';
reKsub = linspace(-1,1,Ksub+1)';
for ii = 1:25
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)];
    Ve = bspline_basismatrix(NB+1, t, reKsub);
    VX = Ve*re;
    VX = VX(:)';    
    res(ii) = norm(VX-VXopt);
end

figure(1)
vv = kron(speye(KK),Vp)*R*W(:,id);
vv = vv/max(abs(vv));
plot(map(rp),vv,'linewidth',2);hold on

plot(VX,VX*0,'x','markersize',8)
plot(VXopt,VXopt*0,'o','markersize',8)

figure(2)
semilogy(res,'o--','linewidth',2)
hold on

