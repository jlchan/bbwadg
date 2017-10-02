clear
NN = 6;
KK = 128;

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
Ksub = 32; % number of splines
N = NB+Ksub-1;

Bq = kron(speye(KK),Vq)*R;
%DBr = kron(spdiag(2./diff(VX)),Dr)^NB;
DBr = kron(KK^NB*speye(KK),Dr^NB);
Brq = kron(speye(KK),Vq)*DBr*R;

M = Bq'*diag(wBq)*Bq;
K = Brq'*diag(wBq)*Brq;

% apply BCs
mm = 1; %mean(M(:));
kk = 1; %mean(K(:));
for i = 1:NB
    id = i;
    M(id,:) = 0; M(:,id) = 0;  M(id,id) = mm;
    K(id,:) = 0; K(:,id) = 0;  K(id,id) = kk;
    id = size(M,2)-i+1;
    M(id,:) = 0; M(:,id) = 0;  M(id,id) = mm;
    K(id,:) = 0; K(:,id) = 0;  K(id,id) = kk;
end
M = M(NB+1:end-NB,NB+1:end-NB);
K = K(NB+1:end-NB,NB+1:end-NB);

[W D] = eig(M,K);
[lam p] = sort(diag(D),'descend');
W = W(:,p);

VX = linspace(-1,1,Ksub+1);
t0 = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
t = t0(:);
re = linspace(-1,1,N+1)';
sk = 1;
tdiff = 1;
while tdiff > 1e-5
    Ve = bspline_basismatrix(NB+1, t, t0);
    tdiff = norm(t-Ve*re);
    t = Ve*re;
    sk = sk + 1;
end
niter = sk
VX = t(NB+1:end-NB)';

%id = 2*NB+Ksub; % skip over BCs
id = Ksub; % skip over BCs
w = [zeros(NB,1);W(:,id);zeros(NB,1)];
    

computeRoots = 1;
if computeRoots
    % compute optimal knots as roots
    VXK = linspace(-1,1,KK+1);
    tK = [VXK(1)*ones(1,NN) VXK VXK(end)*ones(1,NN)]; % open knot vec
    VXopt = zeros(size(VX));
    VXopt(1) = -1;
    VXopt(end) = 1;
    nroots = floor((length(VX)-2)/2)+1;
    for ii = 2:nroots
        disp(sprintf('computing root %d out of %d\n',ii,nroots));
        VXopt(ii) = fzero(@(r) bspline_basismatrix(NN+1, tK, r)*w, VX(ii));
    end
    VXopt(end-nroots+1:end) = -fliplr(VXopt(1:nroots)); % use symmetry        
end

%% plotting and errors with smoothed knots
% plot optimal knots
if 0
    vv = kron(speye(KK),Vp)*R*W(:,id);
    vv = vv/max(abs(vv));
    plot(map(rp),-vv,'b-','linewidth',2);hold on
    plot(VXopt,VXopt*0,'ro','markersize',8,'linewidth',2)
    grid on
    set(gca,'fontsize',14)
    print(gcf,'-dpng','../fnwidth_N2_K8.png')
    % print(gcf,'-dpng','../fnwidth_N3_K8.png')
    % print(gcf,'-dpng','../fnwidth_N4_K8.png')
    return
end

if computeRoots
    % compute errors
    VX = linspace(-1,1,Ksub+1);
    t0 = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]'; % open knot vec
    t = t0;
    re = linspace(-1,1,N+1)';
    for ii = 1:niter
        Ve = bspline_basismatrix(NB+1, t, t0);
        t = Ve*re;
        
        topt = [VXopt(1)*ones(1,NB) VXopt VXopt(end)*ones(1,NB)];
        
        res(ii) = norm(t(:)-topt(:));
        
        if 0
            rp = linspace(-1,1,250)';
            Vp = bspline_basismatrix(NB+1, t, rp);
            Vpopt = bspline_basismatrix(NB+1, topt, rp);
            
            clf
            plot(rp,Vp)
            hold on
            plot(rp,Vpopt,'--')
            %         plot(VX,VX*0,'o','markersize',8,'linewidth',2)
            %         plot(VXopt,VXopt*0,'x','markersize',8,'linewidth',2)
            plot(t0,Ve,'o','markersize',8,'linewidth',2)
            
            plot(t0,(1+t)/2,'x--','markersize',8)
            plot(t0,(1+topt)/2,'d--','markersize',8)
            pause
        end
        
    end
end

figure(1)
vv = kron(speye(KK),Vp)*R*w;
vv = vv/max(abs(vv));
plot(map(rp),vv,'b-','linewidth',2);hold on

if computeRoots
    plot(topt,topt*0,'ro','markersize',10,'linewidth',2)
    %     figure(2)
    %     semilogy(res,'o--','linewidth',2)
    %     hold on
    %     title('Error between opt and smoothed knots with iteration')
    %     ylabel('Smoothing iteration','fontsize',14)
end

plot(t,t*0,'kx','markersize',10,'linewidth',2)
grid on
set(gca,'fontsize',14)
axis([-1 1 -1.5 1.5])
legend('Eigenfunction','Optimal knots','Smoothed knots','Location','NorthWest')

%  print(gcf,'-dpng','../smoothKnot.png')

clear r ropt
for i = 1:NB+Ksub
    r(i,1) = mean(t((i+1):(i+NB)));
    ropt(i,1) = mean(topt((i+1):(i+NB)));
end

max(abs(t(:)-topt(:)))
max(abs(r-ropt))