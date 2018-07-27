% clear
N = 5;
[r1D w1D] = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);

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

a = .0;
dx = sin(pi*(1+r)/2).*sin(pi*(1+s)/2).*sin(pi*(1+t)/2);
dy = 0;
dz = 0;
x = r + a*dx;
y = s + a*dy;
z = t + a*dz;

V1D = Vandermonde1D(N,r1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
I = eye(length(r1D));
Dr = kron(kron(I,D1D),I);
Ds = kron(kron(I,I),D1D);
Dt = kron(kron(D1D,I),I);

rxJ = Dt*((Ds*y).*z) - Ds*((Dt*y).*z); % this is the problematic one
sxJ = Dr*((Dt*y).*z) - Dt*((Dr*y).*z);
txJ = Ds*((Dr*y).*z) - Dr*((Ds*y).*z);

ryJ = -(Dt*((Ds*x).*z) - Ds*((Dt*x).*z));
syJ = -(Dr*((Dt*x).*z) - Dt*((Dr*x).*z));
tyJ = -(Ds*((Dr*x).*z) - Dr*((Ds*x).*z));

rzJ = -(Dt*((Ds*y).*x) - Ds*((Dt*y).*x));
szJ = -(Dr*((Dt*y).*x) - Dt*((Dr*y).*x));
tzJ = -(Ds*((Dr*y).*x) - Dr*((Ds*y).*x));

fprintf('GCL for elem  = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))

% loadVf; % only for N=2


rf = Vf*r;
sf = Vf*s;
tf = Vf*t;
rxJf = Vf*rxJ; sxJf = Vf*sxJ; txJf = Vf*txJ;
ryJf = Vf*ryJ; syJf = Vf*syJ; tyJf = Vf*tyJ;
rzJf = Vf*rzJ; szJf = Vf*szJ; tzJf = Vf*tzJ;

% nxJ = rxJ(fids).*nrJ + sxJ(fids).*nsJ + txJ(fids).*ntJ;
% nyJ = ryJ(fids).*nrJ + syJ(fids).*nsJ + tyJ(fids).*ntJ;
% nzJ = rzJ(fids).*nrJ + szJ(fids).*nsJ + tzJ(fids).*ntJ;

[wr ws] = meshgrid(w1D);
wf = wr(:).*ws(:);
wf = [wf;wf;wf;wf;wf;wf];
% wf = wf.^0;
[wr ws wt] = meshgrid(w1D);
wq = wr(:).*ws(:).*wt(:);
W = diag(wq);

Qr = W*Dr;
Qs = W*Ds;
Qt = W*Dt;

nrJ = zeros(size(Vf,1),1);
nsJ = zeros(size(Vf,1),1);
ntJ = zeros(size(Vf,1),1);
fids = 1:(N+1)^2;
nrJ(fids) = -1; fids = fids + (N+1)^2;
nrJ(fids) = 1; fids = fids + (N+1)^2;
nsJ(fids) = -1; fids = fids + (N+1)^2;
nsJ(fids) = 1; fids = fids + (N+1)^2;
ntJ(fids) = -1; fids = fids + (N+1)^2;
ntJ(fids) = 1;

zzr = 0*diag(nrJ.*wf);
zzs = 0*diag(nsJ.*wf);
zzt = 0*diag(ntJ.*wf);
QNr = [(Qr-Qr') Vf'*diag(nrJ.*wf);
    -diag(nrJ.*wf)*Vf zzr];

QNs = [(Qs-Qs') Vf'*diag(nsJ.*wf);
    -diag(nsJ.*wf)*Vf zzs];

QNt = [(Qt-Qt') Vf'*diag(ntJ.*wf);
    -diag(ntJ.*wf)*Vf zzt];

DNr = diag(1./[wq;wf])*QNr;
DNs = diag(1./[wq;wf])*QNs;
DNt = diag(1./[wq;wf])*QNt;
DNs(abs(DNs)<1e-8) = 0;

norm(DNr(1:(N+1)^3,:)*[rxJ;rxJf] + DNs(1:(N+1)^3,:)*[sxJ;sxJf] + DNt(1:(N+1)^3,:)*[txJ;txJf] ,'fro')
norm(DNr(1:(N+1)^3,:)*[ryJ;ryJf] + DNs(1:(N+1)^3,:)*[syJ;syJf] + DNt(1:(N+1)^3,:)*[tyJ;tyJf] ,'fro')
norm(DNr(1:(N+1)^3,:)*[rzJ;rzJf] + DNs(1:(N+1)^3,:)*[szJ;szJf] + DNt(1:(N+1)^3,:)*[tzJ;tzJf] ,'fro')
% [rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ;rxJf sxJf txJf ryJf syJf tyJf rzJf szJf tzJf]

plot3(r,s,t,'o');
hold on
text(r+.05,s,t,num2str((1:length(r(:)))'))
xlabel('r')
ylabel('s')

[sxJi sxJj] = meshgrid([sxJ;sxJf]);
sxJavg = .5*(sxJi + sxJj);

[syJi syJj] = meshgrid([syJ;syJf]);
syJavg = .5*(syJi + syJj);

[szJi szJj] = meshgrid([szJ;szJf]);
szJavg = .5*(szJi + szJj);

% sum(DNs(1:(N+1)^3,(N+1)^3+1:end).*syJavg(1:(N+1)^3,(N+1)^3+1:end),2)
% [sum(DNs(1:(N+1)^3,(N+1)^3+1:end).*sxJavg(1:(N+1)^3,(N+1)^3+1:end),2),...
%     sum(DNs(1:(N+1)^3,(N+1)^3+1:end).*syJavg(1:(N+1)^3,(N+1)^3+1:end),2),...
%     sum(DNs(1:(N+1)^3,(N+1)^3+1:end).*szJavg(1:(N+1)^3,(N+1)^3+1:end),2)]



