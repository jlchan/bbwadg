clear
N = 4;
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

NODETOL = 1e-8;

[r2 s2] = meshgrid(r1D);
r2 = r2(:);
s2 = s2(:);
e = ones(size(r2));

rf = [-e; e; r2; r2; r2; r2];
sf = [r2; r2; -e; e; s2; s2];
tf = [s2; s2; s2; s2; -e; e];

Vf = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,r,s,t);
NODETOL = 1e-8;

Kvec = [2 4 8 16 32 64];
sk = 1;
for K1D = Kvec 
    h = 2/K1D;
    x = h*r;
    y = h*s;
    z = h*t;
    
    a = 1;
    d = exp(x+y+z+1/3);
    x = x + a*d;
    y = y + a*d;
    z = z + a*d;
    
    V1D = Vandermonde1D(N,r1D);
    V = kron(kron(V1D,V1D),V1D);
    D1D = GradVandermonde1D(N,r1D)/V1D;
    I = eye(length(r1D));
    Dr = kron(kron(I,D1D),I);
    Ds = kron(kron(I,I),D1D);
    Dt = kron(kron(D1D,I),I);        
        
    % filtering
    e = ones(N+1,1); e(end) = 0; % reduce degree by 1
    VNm1 = Vandermonde1D(N-1,JacobiGL(0,0,N-1));
    F1D = (Vandermonde1D(N-1,JacobiGL(0,0,N))/VNm1)*(Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D);
%     F1D = eye(N+1); % no filter
    
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
    Vq1D = Vandermonde1D(N,rq1D)/V1D;
    Vq = kron(kron(Vq1D,Vq1D),Vq1D);
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
    err(sk) = sqrt(abs(K*sum(sum(wJq.*((rxJex-Vq*rxJ).^2+(sxJex-Vq*sxJ).^2+(txJex-Vq*txJ).^2+...
        (ryJex-Vq*ryJ).^2+(syJex-Vq*syJ).^2+(tyJex-Vq*tyJ).^2+...
        (rzJex-Vq*rzJ).^2+(szJex-Vq*szJ).^2+(tzJex-Vq*tzJ))))));
    sk = sk + 1;
end

hvec = 2./Kvec;
loglog(hvec,err,'o--')
hold on
% r = N+2;
% scale = .025*err(1)/(hvec(1).^r);
% loglog(hvec,scale*hvec.^r,'k--')

