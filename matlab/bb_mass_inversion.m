% check bernstein mass matrix inversion

%% 2D test
clear

for N = 1:15    
    [rq sq wq] = Cubature2D(2*N+2);
    Vq = bern_basis_tri(N,rq,sq);
    M = Vq'*diag(wq)*Vq;
    Np = size(M,2);
    
    [r s] = Nodes2D(N); [r s] = xytors(r,s);
    E = {};
    for i = 0:N
        Ei = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
        Ei(abs(Ei)<1e-8) = 0;
        E{i+1} = Ei;
    end
    
    d = 2;  Dhat = sum(wq);
    lam = zeros(N+1,1);
    for i = 0:N
        denom = factorial(N+i+d)*factorial(N-i);
        lam(i+1) = Dhat*factorial(N)^2*factorial(d) / denom;
    end    
    for j = 0:N
        ilamj = zeros(N-j+1,1);
        for i = 0:N-j
            Nj = N-j;
            denom = factorial(Nj+i+d)*factorial(Nj-i);
            ilamj(i+1) = Dhat*factorial(Nj)^2*factorial(d) / denom;
        end
        ilam{j+1} = ilamj;
        L(1:length(ilamj),j+1) = ilamj./lam(1:length(ilamj));
    end
    
    b = 1./lam;
    c = L\b;
    
    P = 0;
    for i = 0:N
        P = P + c(i+1)*E{i+1}*E{i+1}';
    end
    
    invM = inv(M);
    norm(P-invM,'fro')/norm(invM,'fro')
    
    % x = ones(Np,1);
    % b = 4/3 / Np * ones(Np,1);
    x = randn(Np,1);
    b = M*x;
    norm(P*M - eye(Np),'fro')
    
    % test approx on exp(x+y+z)
    f = @(r,s) exp((r+s)/3);
    VDM = Vandermonde2D(N,rq,sq);
    err = f(rq,sq) - VDM*(VDM'*(wq.*f(rq,sq)));
    err1(N) = sqrt(sum(err(:).^2));
    
    b = Vq'*(wq.*f(rq,sq)); 
    Pb = 0; 
    for i = 0:N
        Pb = Pb + c(i+1)*E{i+1}*E{i+1}'*b; 
    end 
%     for i = 0:N
%         bi{i+1} = E{i+1}'*b;
%     end
    err = f(rq,sq) - Vq*(Pb);
    err2(N) = sqrt(sum(err(:).^2));
    
    b = Vq'*(wq.*f(rq,sq));     
    err = f(rq,sq) - Vq*(M\b);
    err3(N) = sqrt(sum(err(:).^2));
end
semilogy(err1,'o--')
hold on
semilogy(err2,'s--')
semilogy(err3,'x--')

%% 3D test

clear

for N = 1:9
    [rq sq tq wq] = tet_cubature(2*N+2);
    Vq = bern_basis_tet(N,rq,sq,tq);
    M = Vq'*diag(wq)*Vq;
    Np = size(M,2);
    
    [r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
    E = {};
    for i = 0:N
        Ei = bern_basis_tet(N,r,s,t)\bern_basis_tet(N-i,r,s,t);
        Ei(abs(Ei)<1e-8) = 0;
        E{i+1} = Ei;
    end
    
    [r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
    Ei = {};
    for i = 1:N
        Ee = bern_basis_tet(N-i+1,r,s,t)\bern_basis_tet(N-i,r,s,t);
        Ee(abs(Ee)<1e-8) = 0;
        Ei{i} = Ee;
    end
    
    d = 3; Dhat = sum(wq);
    lam = zeros(N+1,1);
    for i = 0:N        
        denom = factorial(N+i+d)*factorial(N-i);
        lam(i+1) = Dhat*factorial(N)^2*factorial(d) / denom;
    end
    
    for j = 0:N
        ilamj = zeros(N-j+1,1);
        for i = 0:N-j
            d = 3;
            Nj = N-j;
            denom = factorial(Nj+i+d)*factorial(Nj-i);
            ilamj(i+1) = Dhat*factorial(Nj)^2*factorial(d) / denom;
        end
        ilam{j+1} = ilamj;
        L(1:length(ilamj),j+1) = ilamj./lam(1:length(ilamj));
    end
    
    b = 1./lam;
    c = L\b;
    
    P = 0;
    for i = 0:N
        P = P + c(i+1)*E{i+1}*E{i+1}';
    end
    
    invM = inv(M);
    norm(P-invM,'fro')/norm(invM,'fro')
    
    % x = ones(Np,1);
    % b = 4/3 / Np * ones(Np,1);
    x = randn(Np,1);
    b = M*x;
    norm(P*M - eye(Np),'fro')
    
    % test approx on exp(x+y+z)
    f = @(r,s,t) exp((r+s+t)/3);
    VDM = Vandermonde3D(N,rq,sq,tq);
    err = f(rq,sq,tq) - VDM*(VDM'*(wq.*f(rq,sq,tq)));
    err1(N) = sqrt(sum(err(:).^2));
    
    b = Vq'*(wq.*f(rq,sq,tq));
    
%     Pb = 0;
%     for i = 0:N
%         Pb = Pb + c(i+1)*E{i+1}*E{i+1}'*b;
%     end       
    bi{1} = Ei{1}'*b;
    for i = 2:N
        bi{i} = Ei{i}'*bi{i-1};
    end    
    bn = c(N+1)*bi{N};
    for i = N:-1:2
        bn = c(i)*bi{i-1} + Ei{i}*bn;
    end    
    Pb = c(1)*b + Ei{1}*bn;
    
    err = f(rq,sq,tq) - Vq*(Pb);
    err2(N) = sqrt(sum(err(:).^2));
    
    b = Vq'*(wq.*f(rq,sq,tq));     
    err = f(rq,sq,tq) - Vq*(M\b);
    err3(N) = sqrt(sum(err(:).^2));
end
% [err1(N),err2(N),err3(N)]
% return
semilogy(err1,'o--')
hold on
semilogy(err2,'s--')
semilogy(err3,'x--')
