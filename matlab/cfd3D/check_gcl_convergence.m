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

if 0
    [x,y,~] = meshgrid(linspace(-1,1,10));
    
    a = .25;
    dx = cos(x).*sin(y).*sin(z);
    dy = sin(x).*cos(y).*sin(z);
    dz = sin(x).*sin(y).*cos(z);
    x = x + a*dx;
    y = y + a*dy;
    z = z + a*dz;
    
    surf(x,y,z)
end

if N>1
    VNm1 = Vandermonde1D(N-1,JacobiGL(0,0,N-1));
    F1D = (Vandermonde1D(N-1,JacobiGL(0,0,N))/VNm1)*(Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D);
else
    F1D = .5*ones(2);
end
Filter_r = sparse(kron(kron(I,F1D),I));
Filter_s = sparse(kron(kron(I,I),F1D));
Filter_t = sparse(kron(kron(F1D,I),I));

Kvec = 2.^(1:4);
sk = 1;
err1 = zeros(size(Kvec)); 
err2 = zeros(size(Kvec));
for K1D = Kvec
    for i = 1:K1D
        for j = 1:K1D
            for k = 1:K1D
                h = 2/K1D;
                x = h*(r + 2*i);
                y = h*(s + 2*j);
                z = h*(t + 2*k);
                
                a = .25;
                dx = cos(x).*sin(y).*sin(z);
                dy = sin(x).*cos(y).*sin(z);
                dz = sin(x).*sin(y).*cos(z);
                x = x + a*dx;
                y = y + a*dy;
                z = z + a*dz;                                
                
                Fr = (Dr*y).*z;
                Fs = (Ds*y).*z;
                Ft = (Dt*y).*z;                
                Fr = Filter_r*Fr;
                Fs = Filter_s*Fs;
                Ft = Filter_t*Ft;
                rxJ = Dt*(Fs) - Ds*(Ft);
                sxJ = Dr*(Ft) - Dt*(Fr);
                txJ = Ds*(Fr) - Dr*(Fs);
                
                Fr = (Dr*x).*z;
                Fs = (Ds*x).*z;
                Ft = (Dt*x).*z;
                Fr = Filter_r*Fr;
                Fs = Filter_s*Fs;
                Ft = Filter_t*Ft;
                ryJ = -(Dt*(Fs) - Ds*(Ft));
                syJ = -(Dr*(Ft) - Dt*(Fr));
                tyJ = -(Ds*(Fr) - Dr*(Fs));
                
                Fr = (Dr*y).*x;
                Fs = (Ds*y).*x;
                Ft = (Dt*y).*x;
                Fr = Filter_r*Fr;
                Fs = Filter_s*Fs;
                Ft = Filter_t*Ft;
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
                wJq = wq.*J;
                
                rxJex =  (ys.*zt - zs.*yt); ryJex = -(xs.*zt - zs.*xt); rzJex = (xs.*yt - ys.*xt);
                sxJex = -(yr.*zt - zr.*yt); syJex =  (xr.*zt - zr.*xt); szJex = -(xr.*yt - yr.*xt);
                txJex =  (yr.*zs - zr.*ys); tyJex = -(xr.*zs - zr.*xs); tzJex = (xr.*ys - yr.*xs);                
                
                Gnorm = 1;
                %     Gnorm = sqrt(K*abs(sum(sum(wJq.*(rxJex.^2+sxJex.^2+txJex.^2 + ...
                %         ryJex.^2+syJex.^2+tyJex.^2 + ...
                %         rzJex.^2+szJex.^2+tzJex.^2)))));
                err1(sk) = err1(sk) + abs(sum(sum(wJq.*(...
                    (rxJex-Vq*rxJ).^2+(sxJex-Vq*sxJ).^2+(txJex-Vq*txJ).^2+...
                    (ryJex-Vq*ryJ).^2+(syJex-Vq*syJ).^2+(tyJex-Vq*tyJ).^2+...
                    (rzJex-Vq*rzJ).^2+(szJex-Vq*szJ).^2+(tzJex-Vq*tzJ).^2)...
                    )))/Gnorm;
                
                %     gclerr = norm(Dr*rxJ + Ds*sxJ + Dt*txJ) + norm(Dr*ryJ + Ds*syJ + Dt*tyJ) + norm(Dr*rzJ + Ds*szJ + Dt*tzJ);
                %     fprintf('GCL error Hdiv = %g\n',gclerr)
                
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
                err2(sk) = err2(sk) + abs(sum(sum(wJq.*(...
                    (rxJex-Vq*rxJ).^2+(sxJex-Vq*sxJ).^2+(txJex-Vq*txJ).^2+...
                    (ryJex-Vq*ryJ).^2+(syJex-Vq*syJ).^2+(tyJex-Vq*tyJ).^2+...
                    (rzJex-Vq*rzJ).^2+(szJex-Vq*szJ).^2+(tzJex-Vq*tzJ).^2)...
                    )))/Gnorm;
                
                %     gclerr = norm(Dr*rxJ + Ds*sxJ + Dt*txJ) + norm(Dr*ryJ + Ds*syJ + Dt*tyJ) + norm(Dr*rzJ + Ds*szJ + Dt*tzJ);
                %     fprintf('GCL error 2 = %g\n',gclerr)
                %     fprintf('=====================================\n')
            end
        end
    end    
    sk = sk + 1;
end
err1 = sqrt(err1);
err2 = sqrt(err2);

hvec = 2./Kvec;
C1 = [ones(size(hvec(end-2:end)))' log(hvec(end-2:end))']\log(err1(end-2:end)');
C2 = [ones(size(hvec(end-2:end)))' log(hvec(end-2:end))']\log(err2(end-2:end)');

print_pgf_coordinates(hvec,err1)
print_pgf_coordinates(hvec,err2)
[C1(2) C2(2)]

figure(1)
loglog(hvec,err1,'o-','linewidth',2)
hold on
loglog(hvec,err2,'x--','linewidth',2)
legend('Hdiv','Standard')
title(sprintf('Computed rates = %f, %f\n',C1(2),C2(2)))
