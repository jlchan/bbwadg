clear
clc

sk = 1;
NB = 3;
for K = [8 16 32 64]
    clear wq MM
    smoothKnots = 0;
    rp = linspace(-1,1,125);
    [Vp M Dr R rBq wBq Bq Brq VX] = bsplineVDM(NB,K,rp,smoothKnots);    
    Pq = M\(Bq'*diag(wBq));
        
    knots = VX;
    midpts = VX(1:end-1)+.5*diff(VX);
    
    rq0 = VX(1) + (VX(2)-VX(1))*linspace(0,1,NB+1);
    rqK = VX(end-1) + (VX(end)-VX(end-1))*linspace(0,1,NB+1);
    rq = sort(uniquetol([rq0 rqK knots midpts]));
    rq = rq(:);
    
    [Vq M] = bsplineVDM(NB,K,rq,smoothKnots);
            
    % construct weighted quadrature: M = Vq'*(W*Vq) 
    WVq = pinv(Vq')*M;         
    
    Pq2 = M\(WVq'); % same as pinv(Vq)!    
%     pinv(WVq) = M\Vq';
    
    f = @(x) exp(sin(pi*x));
    w = @(x) 1 + exp(x);        
    
    wproj = (Bq'*diag(w(rBq).*wBq)*Bq)\(Bq'*diag(wBq))*f(rBq);
    wadg = Pq*(diag(1./w(rBq))*Bq*Pq*f(rBq));
    b = WVq'*f(rq);
    
    approxMwinv = inv(M) * WVq'*diag(1./w(rq))*Vq * inv(M);
%     approxMwinv = pinv(Vq)*diag(1./w(rq))* pinv(WVq'); 
    wadg_wq = approxMwinv * b;       
    
    L2err_wproj(sk) = sqrt((wproj-wadg)'*M*(wproj-wadg));
    L2err_wadg(sk) = sqrt((wproj-wadg_wq)'*M*(wproj-wadg_wq));        
    
    L2errp(sk) = sqrt(wBq'*(f(rBq)-Bq*Pq*f(rBq)).^2);
    L2errw(sk) = sqrt(wBq'*(f(rBq)-Bq*Pq2*f(rq)).^2);
    
    h(sk) = 1/K;
    sk = sk + 1;    
end
loglog(h,L2errp,'o--')
hold on
loglog(h,L2errw,'x--')

% loglog(h,L2err_wproj,'^-')
% hold on
% loglog(h,L2err_wadg,'s-')

loglog(h,h.^2,'--')

legend('L2 projection error','Weighted quadrature projection error','h^2')
