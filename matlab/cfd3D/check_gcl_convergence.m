clear

N = 3;
[r1D w1D] = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);
[wr ws wt] = meshgrid(w1D);
w = wr(:).*ws(:).*wt(:);

[rq1D wq1D] = JacobiGQ(0,0,N+2);
[rq sq tq] = meshgrid(r1D);
rq = rq(:); sq = sq(:); tq = tq(:);
[wr ws wt] = meshgrid(wq1D);
wq = wr(:).*ws(:).*wt(:);

V1D = Vandermonde1D(N,r1D);
V = kron(kron(V1D,V1D),V1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
I = eye(length(r1D));
Dr = kron(kron(I,D1D),I);
Ds = kron(kron(I,I),D1D);
Dt = kron(kron(D1D,I),I);
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vq = kron(kron(Vq1D,Vq1D),Vq1D);

[r2 s2] = meshgrid(r1D);
r2 = r2(:); s2 = s2(:);
e = ones(size(r2));

rf = [-e; e; r2; r2; r2; r2];
sf = [r2; r2; -e; e; s2; s2];
tf = [s2; s2; s2; s2; -e; e];

Vf = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,r,s,t);
NODETOL = 1e-8;

Kvec = 2.^(2:7);
sk = 1;
for K1D = Kvec 
    h = 2/K1D;
    x = h*r;
    y = h*s;
    z = h*t;
    
%     a = .125;
%     dx = cos(pi/2*x).*sin(pi*y).*sin(pi*z);
%     x = x + a*dx;
%     dy = sin(pi*x).*cos(pi/2*y).*sin(pi*z);
%     y = y + a*dy;
%     dz = sin(pi*x).*sin(pi*y).*cos(pi/2*z);
%     z = z + a*dz;

a = 0.25;
% d = (1+x).*(1-x).*(1-y).*(1+y).*(1+z).*(1-z);
    dx = (1-y).*(1+y).*(1+z).*(1-z);
    dy = (1+x).*(1-x).*(1+z).*(1-z);
    dz = (1+x).*(1-x).*(1-y).*(1+y);
% dx = d;     dy = d;     dz = d;
x = x + a*dx;
y = y + a*dy;
z = z + a*dz;
    
    % a = .1;
    % dx = exp(x+y+z);
    % dy = exp(x+y+z);
    % dz = exp(x+y+z);
    % x = x + a*dx;
    % y = y + a*dy;
    % z = z + a*dz;

%     plot3(x,y,z,'o'); return        
        
    e = ones(N+1,1); e(end) = 0; % reduce degree by 1
    VNm1 = Vandermonde1D(N-1,JacobiGL(0,0,N-1));
    F1D = (Vandermonde1D(N-1,JacobiGL(0,0,N))/VNm1)*(Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D);    
    
    Fr = (Dr*y).*z;
    Fs = (Ds*y).*z;
    Ft = (Dt*y).*z;
    Fr = kron(kron(I,F1D),I)*Fr;
    Fs = kron(kron(I,I),F1D)*Fs;
    Ft = kron(kron(F1D,I),I)*Ft;
    rxJ = Dt*(Fs) - Ds*(Ft);
    sxJ = Dr*(Ft) - Dt*(Fr);
    txJ = Ds*(Fr) - Dr*(Fs);
    
    Fr = (Dr*x).*z;
    Fs = (Ds*x).*z;
    Ft = (Dt*x).*z;
    Fr = kron(kron(I,F1D),I)*Fr;
    Fs = kron(kron(I,I),F1D)*Fs;
    Ft = kron(kron(F1D,I),I)*Ft;
    ryJ = -(Dt*(Fs) - Ds*(Ft));
    syJ = -(Dr*(Ft) - Dt*(Fr));
    tyJ = -(Ds*(Fr) - Dr*(Fs));
        
    Fr = (Dr*y).*x;
    Fs = (Ds*y).*x;
    Ft = (Dt*y).*x;
    Fr = kron(kron(I,F1D),I)*Fr;
    Fs = kron(kron(I,I),F1D)*Fs;
    Ft = kron(kron(F1D,I),I)*Ft;
    rzJ = -(Dt*(Fs) - Ds*(Ft));
    szJ = -(Dr*(Ft) - Dt*(Fr));
    tzJ = -(Ds*(Fr) - Dr*(Fs));    
    
    % est same error on each elem    
    xr = Vq*Dr*x; xs = Vq*Ds*x; xt = Vq*Dt*x;
    yr = Vq*Dr*y; ys = Vq*Ds*y; yt = Vq*Dt*y;
    zr = Vq*Dr*z; zs = Vq*Ds*z; zt = Vq*Dt*z;
    J = xr.*(ys.*zt-zs.*yt) - ...
        yr.*(xs.*zt-zs.*xt) + ...
        zr.*(xs.*yt-ys.*xt);
    rxJex =  (ys.*zt - zs.*yt); ryJex = -(xs.*zt - zs.*xt); rzJex = (xs.*yt - ys.*xt);
    sxJex = -(yr.*zt - zr.*yt); syJex =  (xr.*zt - zr.*xt); szJex = -(xr.*yt - yr.*xt);
    txJex =  (yr.*zs - zr.*ys); tyJex = -(xr.*zs - zr.*xs); tzJex = (xr.*ys - yr.*xs);
    K = K1D.^3;     
    wJq = wq.*(J);
    
    Gnorm = 1; 
%     Gnorm = sqrt(K*abs(sum(sum(wJq.*(rxJex.^2+sxJex.^2+txJex.^2 + ...
%         ryJex.^2+syJex.^2+tyJex.^2 + ...
%         rzJex.^2+szJex.^2+tzJex.^2)))));
    err1(sk) = sqrt(K*abs(sum(sum(wJq.*(...
        (rxJex-Vq*rxJ).^2+(sxJex-Vq*sxJ).^2+(txJex-Vq*txJ).^2+...
        (ryJex-Vq*ryJ).^2+(syJex-Vq*syJ).^2+(tyJex-Vq*tyJ).^2+...
        (rzJex-Vq*rzJ).^2+(szJex-Vq*szJ).^2+(tzJex-Vq*tzJ).^2)...
        ))))/Gnorm;
    
    gclerr = norm(Dr*rxJ + Ds*sxJ + Dt*txJ) + norm(Dr*ryJ + Ds*syJ + Dt*tyJ) + norm(Dr*rzJ + Ds*szJ + Dt*tzJ);
    fprintf('GCL error Hdiv = %g\n',gclerr)
    
    %==========================================
        
    % no filtering
    Fr = (Dr*y).*z;
    Fs = (Ds*y).*z;
    Ft = (Dt*y).*z;
    rxJ = Dt*(Fs) - Ds*(Ft);
    sxJ = Dr*(Ft) - Dt*(Fr);
    txJ = Ds*(Fr) - Dr*(Fs);
    
    Fr = (Dr*x).*z;
    Fs = (Ds*x).*z;
    Ft = (Dt*x).*z;
    ryJ = -(Dt*(Fs) - Ds*(Ft));
    syJ = -(Dr*(Ft) - Dt*(Fr));
    tyJ = -(Ds*(Fr) - Dr*(Fs));
        
    Fr = (Dr*y).*x;
    Fs = (Ds*y).*x;
    Ft = (Dt*y).*x;
    rzJ = -(Dt*(Fs) - Ds*(Ft));
    szJ = -(Dr*(Ft) - Dt*(Fr));
    tzJ = -(Ds*(Fr) - Dr*(Fs));   
    err2(sk) = sqrt(K*abs(sum(sum(wJq.*(...
        (rxJex-Vq*rxJ).^2+(sxJex-Vq*sxJ).^2+(txJex-Vq*txJ).^2+...
        (ryJex-Vq*ryJ).^2+(syJex-Vq*syJ).^2+(tyJex-Vq*tyJ).^2+...
        (rzJex-Vq*rzJ).^2+(szJex-Vq*szJ).^2+(tzJex-Vq*tzJ).^2)...
        ))))/Gnorm;
    
    gclerr = norm(Dr*rxJ + Ds*sxJ + Dt*txJ) + norm(Dr*ryJ + Ds*syJ + Dt*tyJ) + norm(Dr*rzJ + Ds*szJ + Dt*tzJ);
    fprintf('GCL error 2 = %g\n',gclerr)
    fprintf('=====================================\n')
    sk = sk + 1;
end

hvec = 2./Kvec;
C1 = [ones(size(hvec(end-2:end)))' log(hvec(end-2:end))']\log(err1(end-2:end)');
C2 = [ones(size(hvec(end-2:end)))' log(hvec(end-2:end))']\log(err2(end-2:end)');

loglog(hvec,err1,'o-','linewidth',2)
hold on
loglog(hvec,err2,'x--','linewidth',2)
legend('Hdiv','Standard')
title(sprintf('Computed rates = %f, %f\n',C1(2),C2(2)))
