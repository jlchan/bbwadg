clear

N = 3;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
V = Vandermonde2D(N,r,s);
[Vr Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/V;
Ds = Vs/V;

[rq sq wq] = Cubature2D(2*N);
Vq = Vandermonde2D(N,rq,sq)/V;
Vqref = Vq;

[rq1D wq1D] = JacobiGQ(0,0,N+3);
e = ones(length(rq1D),1);
rf = [rq1D;-rq1D;-e];
sf = [-e;rq1D;-rq1D];
wf = [wq1D;wq1D;wq1D];
nrJ = [0*e;e;-e];
nsJ = [-e;e;0*e];
Vf = Vandermonde2D(N,rf,sf)/V;

%% make physical mapping

k = pi;
a = .1;
x = r + a*sin(k*(.1+r+s));
y = s + a*sin(k*(.1+r+s));

[rx,ry,sx,sy,J] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);

[rxf,ryf,sxf,syf,Jf] = GeometricFactors2D(x,y,Vf*Dr,Vf*Ds);
rxJ = rxf.*Jf; sxJ = sxf.*Jf;
ryJ = ryf.*Jf; syJ = syf.*Jf;
nxJ = rxJ.*nrJ + sxJ.*nsJ;
nyJ = ryJ.*nrJ + syJ.*nsJ;

sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ;
ny = nyJ./sJ;

xq = Vq*x; yq = Vq*y;
xf = Vf*x; yf = Vf*y;

plot(xq,yq,'o')
hold on
plot(xf,yf,'x')
quiver(xf,yf,nx,ny)
% return

Vphys = zeros(length(xq),(N+1)*(N+2)/2);
Vfphys = zeros(length(xf),(N+1)*(N+2)/2);
Vxphys = zeros(length(xq),(N+1)*(N+2)/2);
Vyphys = zeros(length(xq),(N+1)*(N+2)/2);
sk = 1;
for i = 0:N    
    for j = 0:N-i
        Vphys(:,sk) = xq.^i .* yq.^j;
        Vfphys(:,sk) = xf.^i .* yf.^j;
        if i > 0
            Vxphys(:,sk) = i*xq.^(i-1).*yq.^j;
        end
        if j > 0
            Vyphys(:,sk) = xq.^i.*j.*yq.^(j-1);
        end
        sk = sk + 1;
    end
end

% Q: modes->qpts, Vphys: monomials->qpts
% R = monomials->modes, inv(R): modes->monomials
[Q R] = qr(Vphys,0);  % Q*R = Vphys 
Vq = Q;
Vf = Vfphys/R;
M = Vq'*diag(wq.*J)*Vq;

Pq = M\(Vq'*diag(wq.*J)); % Pq: qpts -> modes

E = Vf*Pq; 

Dx = Pq*(Vxphys/R); 
Dy = Pq*(Vyphys/R); 
Qx = diag(wq.*J)*Vq*Dx*Pq;
Qy = diag(wq.*J)*Vq*Dy*Pq;

norm((diag(1./(wq.*J))*Qx)*xq.^N-N*xq.^(N-1))

Vqref*(Vqref\(xq+yq))-(xq+yq)
return

e = ones(size(Qx,2),1);
% norm(Qx'*e - E'*diag(nxJ)*wf,'fro')

% solve min .5*||wf - wfnew||^2 st exactness constraints
Ax = Vq'*E'*diag(nxJ);
Ay = Vq'*E'*diag(nyJ);
A = [Ax;Ay];

[UA SA VA] = svd(A,0);
nnz(diag(SA)>1e-14) 

C = [eye(length(wf)) A';
    A zeros(size(A,1))];
b = [wf; [Vq'*Qx'*e;Vq'*Qy'*e]];
ww = C\b;
wfnew = ww(1:length(wf));
% bxy = [Qx'*e;Qy'*e] - A*wf;
% wfnew = A'*pinv(A*A')*bxy + wf;
norm(Qx'*e - E'*diag(nxJ)*wfnew,'fro')
norm(Qy'*e - E'*diag(nyJ)*wfnew,'fro')
norm(wf - wfnew)

% (Vfphys/R)'*diag(nx)*(wf.*sJ) - (Vxphys/R)'*Vq'*(wq.*J)

% vv = J;
% color_line3(xq,yq,vv,vv,'.')


