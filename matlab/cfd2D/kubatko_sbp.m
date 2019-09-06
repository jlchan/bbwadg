load ~/Downloads/QuadratureRules.mat

N = 6;
quadpts = Q_GaussLegendre{N};
[r1D,w1D] = JacobiGQ(0,0,N);
% quadpts = Q_GaussLobatto{N};
% [r1D,w1D] = JacobiGL(0,0,N+1); % 2*(N+1)-1 = 2*N+1

rs = quadpts.Points;
wq = quadpts.Weights;

rq = rs(:,1);
sq = rs(:,2);

% vol nodes
Vq = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Dr = Vq\Vr; Ds = Vq\Vs;

M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

% face nodes
ef = ones(size(r1D));
rf = [r1D; r1D; -ef];
sf = [-ef; -r1D; r1D];
wf = [w1D; w1D; w1D];
Vf = Vandermonde2D(N,rf,sf);
Ef = Vf*Pq;

E = zeros(length(rf),length(rq));
for i = 1:length(rf)
    fid = find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10);
    E(i,fid) = 1;
end

e = ones(length(rf)/3,1);
nrJ = [0*e; e; -e];
nsJ = [-e; e; 0*e];

Qr = diag(wq)*Vq*Dr*Pq;
Qs = diag(wq)*Vq*Ds*Pq;


VN = [eye(length(wq));E]; 
Br = diag(nrJ.*wf);
Bs = diag(nsJ.*wf);
QNr = [Qr-.5*Ef'*Br*Ef .5*Ef'*Br;
          -.5*Br*Ef .5*Br];
QNs = [Qs-.5*Ef'*Bs*Ef .5*Ef'*Bs;
          -.5*Bs*Ef .5*Bs];
      
Qrsbp = VN'*QNr*VN; 
Qssbp = VN'*QNs*VN; 

% % alternative construction
% QNrskew = [.5*(Qr-Qr') .5*Ef'*Br;
%           -.5*Br*Ef 0*Br];
% QNsskew = [.5*(Qs-Qs') .5*Ef'*Bs;
%           -.5*Bs*Ef 0*Bs];
% BNr = [zeros(length(wq)) zeros(length(wq),length(wf));
%     zeros(length(wf),length(wq)) Br];
% BNs = [zeros(length(wq)) zeros(length(wq),length(wf));
%     zeros(length(wf),length(wq)) Bs];
% Qrsbp = (VN'*(.5*BNr + QNrskew)*VN);
% Qssbp = (VN'*(.5*BNs + QNsskew)*VN);

norm(diag(1./wq)*Qrsbp*((rq+sq).^N) - N*(rq+sq).^(N-1))
norm(Qrsbp+Qrsbp' - E'*diag(wf.*nrJ)*E)

norm(diag(1./wq)*Qssbp*((rq+sq).^N) - N*(rq+sq).^(N-1))
norm(Qssbp+Qssbp' - E'*diag(wf.*nsJ)*E)
% 
% plot(r,s,'x')
% hold on
% plot(rf,sf,'o')
% quiver(rf,sf,nrJ,nsJ)


