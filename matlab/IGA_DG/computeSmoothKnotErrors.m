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

for NB = 2:4 % degree of spline space
    for Ksub = [4 8 16 32 64] % number of splines
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
        while tdiff > 1e-6
            Ve = bspline_basismatrix(NB+1, t, t0);
            tdiff = norm(t-Ve*re);
            t = Ve*re;
            sk = sk + 1;
        end
        niter = sk;
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
%                 disp(sprintf('computing root %d out of %d\n',ii,nroots));
                VXopt(ii) = fzero(@(r) bspline_basismatrix(NN+1, tK, r)*w, VX(ii));
            end
            VXopt(end-nroots+1:end) = -fliplr(VXopt(1:nroots)); % use symmetry
        end
        topt = [VXopt(1)*ones(1,NB) VXopt VXopt(end)*ones(1,NB)];
        
        % compute greville too
        clear r0 rorig ropt
        
        for i = 1:NB+Ksub
            r0(i,1) = mean(t((i+1):(i+NB)));
            rorig(i,1) = mean(t0((i+1):(i+NB)));
            ropt(i,1) = mean(topt((i+1):(i+NB)));
        end
        
        disp(sprintf('N = %d, K = %d: knot error = %f, , original knot error = %f, Greville error = %f, original Greville error = %f\n',...
            NB,Ksub,max(abs(t(:)-topt(:))),max(abs(t0(:)-topt(:))),max(abs(r0-ropt)),max(abs(rorig-ropt))))
        
    end
end
