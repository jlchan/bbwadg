N = 5;
Nq = 3*N;
[rq sq w] = Cubature2D(Nq);

[r1D, w1D] = JacobiGQ(0, 0, Nq);

rfq = [r1D, -r1D, -ones(size(r1D))];
sfq = [-ones(size(r1D)), r1D, -r1D];
wf = [w1D, w1D*sqrt(2), w1D];
rfq = rfq(:);
sfq = sfq(:);
wf = wf(:);
Nfq = length(wf)/3;

Vq = Vandermonde2D(N,rq,sq);
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq; Vsq = Vsq;
Vfq = Vandermonde2D(N,rfq,sfq);

nr = [zeros(size(r1D)) ones(size(r1D)) -ones(size(r1D))];
ns = [-ones(size(r1D)) ones(size(r1D)) zeros(size(r1D))];       
        
sJ = sqrt(nr.^2 + ns.^2);
nr = nr./sJ; ns = ns./sJ;
nr = nr(:); ns = ns(:); 
%sJ = sJ(:);

% plot(rq,sq,'*');hold on
% plot(rfq,sfq,'o');
% hold on;quiver(rfq,sfq,nr,ns)
% text(rfq+.1,sfq,num2str((1:length(rfq))'))
% axis equal
% return

% int(dudr*v) = int(-u*dudr) + int_dK (u * nr * v)
Np = size(Vq,2);
u = (1:Np)'/Np; 
% u = ones(Np,1);
% u = zeros(Np,1); u(1) = 1;
% u = randn(size(Vq,2),1);
norm(Vq'*(w.*(Vrq*u)) - (-Vrq'*(w.*(Vq*u)) + Vfq'*(wf.*nr.*(Vfq*u))),'fro')
