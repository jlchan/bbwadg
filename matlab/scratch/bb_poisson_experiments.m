%% 1D laplacian

N = 20;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
rp = linspace(-1,1,50);
re = linspace(-1,1,N+1)';

[V Vr] = bern_basis_1D(N,r);
Vq = bern_basis_1D(N,rq);
Vp = bern_basis_1D(N,rp);
Dr = V\Vr;
[Vb Vrb] = bern_basis_1D(N,[-1 1]);

% V = Vandermonde1D(N,r);
% Dr = GradVandermonde1D(N,r)/V;
% Vq = Vandermonde1D(N,rq)/V;
% Vp = Vandermonde1D(N,rp)/V;
% Vb = Vandermonde1D(N,[-1 1])/V;
% Vrb = GradVandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;

b = Vq'* (wq.*pi^2.*sin(pi*rq)); 
K = Dr'*M*Dr;

if 1
    K(1,:) = 0; K(:,1) = 0; K(1,1) = 1;
    K(N+1,:) = 0; K(:,N+1) = 0; K(N+1,N+1) = 1;
    b([1 N+1]) = 0;
else % nitsche    
    B = Vrb'*Vb - Vb'*Vrb;
%     e = zeros(N+1,1); e(1) = 1; e(N+1) = 1; e = diag(e);
    %     B = -Vrb'*Vb - Vb'*Vrb + 1e7*Vb'*Vb; % sym nitsche
    K = K + B; % + e*e';
end

invK = inv(K);
rr = JacobiGL(0,0,N+1);
Ep = bern_basis_1D(N+1,rr)\bern_basis_1D(N,rr);
Em = bern_basis_1D(N,rr)\bern_basis_1D(N-1,rr);

plot(rp,Vp*(K\b),'.');
hold on;
plot(rp,sin(pi*rp))
max(abs(sin(pi*rp')-Vp*(K\b)))
return

% invK = invK*diag(1./max(invK,[],1));
% for i = 2:N
%     plot(re,invK(:,i),'o-')
%     axis([-1 1 -1 1])
%     pause
% end
% return

invKp = Vp*invK;
invKp = invKp*diag(1./max(invKp,[],1));
invK = invK*diag(1./max(invK,[],1));
% M = N;
a = 1/(N+1);
% rc = a*cos((2*(0:M))/(2*M)*pi)';
% rc = rc + (1-a)*cos((2*(1:M+1)-1)/(2*M+2)*pi)';
rc = (1-a)*JacobiGL(0,0,N) + a*JacobiGQ(0,0,N);

hold on
plot(rp,invKp)

plot(rc,ones(size(rc)),'o')
[~,ids] = max(invKp(:,2:N));
rpmax = rp(ids);
plot(rpmax,ones(N-1,1),'x')


%% triangle

N = 15;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq wq] = Cubature2D(2*N);

[V Vr Vs] = bern_basis_tri(N,r,s);
Dr = V\Vr; Ds = V\Vs;
Vq = bern_basis_tri(N,rq,sq);

V = Vandermonde2D(N,r,s);
[Vr Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/V; Ds = Vs/V;
Vq = Vandermonde2D(N,rq,sq);

M = Vq'*diag(wq)*Vq;

K = Dr'*M*Dr + Ds'*M*Ds;
bmask = [find(abs(r+1)<1e-8) find(abs(s+1)<1e-8) find(abs(r+s)<1e-8)];
bmask = unique(bmask);
K(bmask,:) = 0;
K(:,bmask) = 0;
K(bmask,bmask) = eye(length(bmask));

iids = setdiff(1:length(r),bmask);
invK = inv(K);
[re se] = EquiNodes2D(N); [re se ] = xytors(re,se);

invK = invK(iids,iids);
K = K(iids,iids);

% subplot(1,2,1)
% plot(re(iids),se(iids),'o');text(re(iids)+.1/N,se(iids),num2str((1:length(iids))'))
% subplot(1,2,2)
% imagesc(invK(iids,iids))

[re se] = EquiNodes2D(N);[re se ] = xytors(re,se);

plot(re,se,'o'); 
idn = 18;
hold on; plot(re(idn),se(idn),'*')
idr = find(abs(Dr(idn,:))>1e-8);
ids = find(abs(Ds(idn,:))>1e-8);
hold on;plot3(re(idr),se(idr),Dr(idn,idr),'^','markersize',20)
hold on;plot3(re(ids),se(ids),Ds(idn,ids),'s','markersize',20)


%% check relation b/w blocks

r1D = JacobiGL(0,0,N);
E = bern_basis_1D(N-3,r1D)\bern_basis_1D(N-4,r1D);
[V1 V1r] = bern_basis_1D(N-3,r1D);
[V2 V2r] = bern_basis_1D(N-4,r1D); 
DE = V1\V2r;

K1 = K(1:N-2,1:N-2);
K2 = K(1:N-2,N-2 + (1:N-3));
K3 = K(1:N-2,N-2 + N-3 + (1:N-4));

[re se] = EquiNodes2D(N-3);[re se ] = xytors(re,se);
plot3(re,se,invK(:,3),'o')

iK1 = invK(1:N-2,1:N-2);
iK2 = invK(1:N-2,N-2 + (1:N-3));
iK3 = invK(1:N-2,N-2 + N-3 + (1:N-4));


KE = iK1\iK2;
for i = 1:size(KE,1)
    for j = 1:size(KE,2)
        if abs(i-j)>2
            KE(i,j) = 0;
        end
    end
end

