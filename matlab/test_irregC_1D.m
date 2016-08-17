clear
Globals1D

N = 3;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,3*N);

avec = [1e-1 1e-4];
b = .1;
adaptQ = 0;

Kvec = [4 8 16 32];
for a = avec
    
    sk = 1;
    
    for K = Kvec
        [Nv, VX, K, EToV] = MeshGen1D(-1,1,K);
        
        StartUp1D
               
        Vq = Vandermonde1D(N,rq)/V;
        xq = Vq*x;
        Jq = Vq*J;        
        
        Vq = Vandermonde1D(N,rq); % use modal for approx
        
        cfun = @(x) 1 + sqrt((x+b).^2 + a);
        ufun = @(x) exp(x);
        
        cq = cfun(xq);
        uq = ufun(xq);
        u1 = [];
        for e = 1:K
            VX1 = VX(EToV(e,1));
            VX2 = VX(EToV(e,2));
            h = VX2-VX1;
            b1 = Vq'*(wq.*uq(:,e));
            if adaptQ
                for i = 0:N
                    for j = 0:N
                        funM = @(r) (JacobiP(r,0,0,i).*JacobiP(r,0,0,j)./cfun((1+r(:))/2*h + VX1))';
                        M1(i+1,j+1) = quadgk(funM,-1,1);
                    end
                end
            else
                M1 = Vq'*diag(wq./cq(:,e))*Vq; 
            end
            
            u1(:,e) = M1\b1;
            
        end
        if adaptQ
            errK = 0;
            for e = 1:K
                VX1 = VX(EToV(e,1));
                VX2 = VX(EToV(e,2));
                h = VX2-VX1;
                funM = @(r) (((Vandermonde1D(N,r)*u1(:,e))./cfun((1+r(:))/2*h + VX1) - ufun((1+r(:))/2*h + VX1)).^2)';                
                errK = errK + J(1,e)*quadgk(funM,-1,1);
            end
            err1(sk) = sqrt(errK);
        else
            ediff = diag(wq)*(Jq.*((Vq*u1)./cq - ufun(xq)).^2);
            err1(sk) = sqrt(sum(ediff(:)));
        end
        sk = sk + 1;
        
    end
    h = .5.^(1:length(err1));
    if adaptQ
        loglog(h,err1,'x-')
    else
        loglog(h,err1,'o-')
    end
    hold on
end


