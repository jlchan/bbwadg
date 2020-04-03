clear
N = 2;

[r s t] = pyr_nodes(N);
[rq sq tq wq] = pyr_cubature(3*N);
[rp sp tp] = pyr_equi_nodes(50);
% tp = linspace(-1,1,100)';
% rp = -1*ones(size(tp));
% sp = -1*ones(size(tp));

[V Vr Vs Vt] = pyr_basis(N,r,s,t);
[Vq Vrq Vsq Vtq] = pyr_basis(N,rq,sq,tq);
Vp = pyr_basis(N,rp,sp,tp)/V;

% [V Vr Vs Vt] = warburton_pyr_basis(N,r,s,t);
% [Vq Vrq Vsq Vtq] = warburton_pyr_basis(N,rq,sq,tq);
% Vp = warburton_pyr_basis(N,rp,sp,tp)/V;

[VB] = bern_pyr(N,r,s,t);

color_line3(rp,sp,tp,Vp(:,1),'.')
view(3)
return
% [r s t] = tet_nodes(N)
% [rq sq tq wq] = tet_cubature(2*N);
% [V Vr Vs Vt] = tet_basis(N,r,s,t);
% [Vq Vrq Vsq Vtq] = tet_basis(N,rq,sq,tq);

% roundabout way of constructing nodal deriv. matrices
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
Drq = Vrq*Pq;
Dsq = Vsq*Pq;
Dtq = Vtq*Pq;
Vqnodal = Vq/V;
Pqnodal = (Vqnodal'*diag(wq)*Vqnodal)\(Vqnodal'*diag(wq));
Dr = Pqnodal*Drq*Vqnodal;
Ds = Pqnodal*Dsq*Vqnodal;
Dt = Pqnodal*Dtq*Vqnodal;

for i = 0:N
    for j = 0:N-i
        for k = 0:N-i-j
           err = norm(Dr*(r.^i.*s.^j.*t.^k) - i*r.^(i-1).*s.^j.*t.^k) + ...
            + norm(Ds*(r.^i.*s.^j.*t.^k) - r.^i.*j.*s.^(j-1).*t.^k) ...
            + norm(Dt*(r.^i.*s.^j.*t.^k) - r.^i.*s.^j.*k.*t.^(k-1));
        end
    end
end

fprintf('\n\n')
fprintf('err in derivs of PN = %g\n',err)
fprintf('derivatives in approx space = %g\n',norm(V*(V\Vr) - Vr,'fro')+norm(V*(V\Vs) - Vs,'fro')+norm(V*(V\Vt) - Vt,'fro'))

if 0
    a = .25;
    x = r;
    y = s;
    z = t;
    x = x + a*cos(x).*sin(y).*sin(z);
    y = y + a*sin(x).*cos(y).*sin(z);
    z = z + a*sin(x).*sin(y).*cos(z);
else
    a = .25;
    [r1 s1 t1] = pyr_nodes(1);
    V1 = pyr_basis(1,r1,s1,t1);
    VN = pyr_basis(1,r,s,t)/V1;       
    x = VN*(r1+a*randn(size(r1)));
    y = VN*(s1+a*randn(size(r1)));
    z = VN*(t1+a*randn(size(r1)));
end

plot3(r,s,t,'o')
hold on
plot3(x,y,z,'o')

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

fprintf('gcl err = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))
fprintf('commuting DrDt err = %g\n',norm(Dr*Dt-Dt*Dr,'fro'))

% norm(Dr,'fro')
% norm(Ds,'fro')
% norm(Dt,'fro')
