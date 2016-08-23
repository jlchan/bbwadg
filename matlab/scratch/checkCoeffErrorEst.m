clear
% check difference in control points vs interpolant
N = 5;
r = JacobiGL(0,0,N);

VB = bern_basis_1D(N,r);

rp = linspace(-1,1,500);
Vp = bern_basis_1D(N,rp);
req = linspace(-1,1,N+1); req = req(:);
[rq w] = JacobiGQ(0,0,5*N);
Vq = bern_basis_1D(N,rq);
% plot(rq,Vq,'-');return
M = Vq'*diag(w)*Vq;

Ve = bern_basis_1D(N,req);
sk = 1;
kvec = .1:.1:10;
for k = kvec
    % fk = sin(k*pi*rq);
    fk = sqrt(rq.^2+k);
    u = M\(Vq'*(w.*fk));
    err(sk) = sqrt(sum(w.*(Vq*u-fk).^2));
    pterr1(sk) = sum(abs(Ve*u-u));
    pterr2(sk) = sqrt(sum((Ve*u-u).^2));
    pterrInf(sk) = max(abs(Ve*u-u));
    sk = sk + 1;
end
semilogy(kvec,err)
hold on;
semilogy(kvec,pterr1,'o--')
semilogy(kvec,pterr2,'--')
semilogy(kvec,pterrInf,'.-')
