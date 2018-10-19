clear
Kvec = [16 32 64 128 256 512 1024];
sk = 1;
for K = Kvec
    
    [H D tL tR x] = HGTpEQ3(K);
    Ivf = [tL tR]';
    B = diag([-1;1]);
    Q = H*D;
    
    norm(Q+Q' - Ivf'*B*Ivf,'fro')
    QN = [Q - .5*Ivf'*B*Ivf .5*Ivf'*B;
        -.5*B*Ivf .5*B];
    
    PN = H\[eye(size(H)) Ivf'];
    
    f = @(x) exp(sin(1+3*pi*x));
    df = @(x) 3*pi*exp(sin(3*pi*x + 1)).*cos(3*pi*x + 1);
    
    %     plot(x,D*f(x),'x--')
    %     hold on;
    %     plot(x,PN*QN*f([x;-1;1]),'o--')
    %     xp = linspace(-1,1,1000)';
    %     plot(xp,df(xp),'-')
    %     return
    
    err1(sk) = max(abs(df(x)-D*f(x)));
    err2(sk) = max(abs(df(x)-PN*QN*f([x;-1;1])));
    sk = sk + 1;
end

h = 1./Kvec;
h = h/h(1);
loglog(h,err1,'o--','linewidth',3)
hold on
loglog(h,err2,'x--','linewidth',3)

loglog(h,err2(2)*h.^3,'k--','linewidth',3)


%%

clear
N = 4;

[r w] = JacobiGQ(0,0,N);
V = Vandermonde1D(N,r);
D = GradVandermonde1D(N,r)/V;
Ivf = Vandermonde1D(N,[-1;1])/V;
H = diag(w);
B = diag([-1;1]);
Q = H*D;

norm(Q+Q' - Ivf'*B*Ivf,'fro')
QN = [Q - .5*Ivf'*B*Ivf .5*Ivf'*B;
    -.5*B*Ivf .5*B];

Ivs = eye(N+1);
Ivs = V*diag([ones(N,1); 0])/V;
PN = H\[Ivs Ivf'];
rN = [r;-1;1];


norm(PN*QN*rN.^N - N*r.^(N-1))


