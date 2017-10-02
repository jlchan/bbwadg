% derivative operators for 

function [Sr Ss St] = pyr_deriv_ops(N)

[c wc] = JacobiGQ(2,0,N);
D = GradVandermonde1D(N,c)/Vandermonde1D(N,c);

a = cell(N,1); b = cell(N,1); w = cell(N,1);
for k = 0:N
    [bk ak wabk] = QCubature2D(k);    
    w{k+1} = wabk;    a{k+1} = ak;    b{k+1} = bk;
end

% a b Vandermonde matrices
VMDab = cell(N,N); Da = cell(N,N); Db = cell(N,N);
for k = 0:N
    VDM = QVandermonde2D(k,a{k+1},b{k+1}); 
    for n = 0:N        
        [Vab Va Vb] = QVandermonde2D(k,a{n+1},b{n+1});        
        % order k nodal basis evaluated at a^n, b^n
        Vab = Vab/VDM; DVa = Va/VDM; DVb = Vb/VDM;
        VDMab{n+1,k+1} = Vab; 
        Da{n+1,k+1} = DVa;
        Db{n+1,k+1} = DVb;
    end
end

% c "mass matrix"
for n = 0:N
    CNn = (N+2)/(2^(2*n+2)*(2*n+3));
    for k = 0:N
        CNk = (N+2)/(2^(2*k+2)*(2*k+3));
        uc = JacobiP(c,2*k+3,0,N-k).*((1-c)/2).^k/sqrt(CNk); 
        vc = JacobiP(c,2*n+3,0,N-n).*((1-c)/2).^n/sqrt(CNn);
        Mc(n+1,k+1) = vc'*((wc.*(2./(1-c))).*uc); % d(ab)/d(rs) factor
        
        % deriv matrix in C
        if k > 0
            duc = GradJacobiP(c,2*k+3,0,N-k).*((1-c)/2).^k + ...
                -.5*k*((1-c)/2).^(k-1).*JacobiP(c,2*k+3,0,N-k);
        else
            duc = GradJacobiP(c,2*k+3,0,N-k);
        end 
        duc = duc/sqrt(CNk);
        S1Dc(n+1,k+1) = vc'*(wc.*duc); % d(ab)/d(rs) factor
    end
end
Mc(abs(Mc)<1e-12) = 0;

Np = (N+1)*(N+2)*(2*N+3)/6;

Sr = zeros(Np);
rsk = 1;
for n = 1:N+1
    for lm = 1:n^2
        csk = 1;
        for k = 1:N+1
            for ij = 1:k^2
                sqw = sqrt(w{n}(lm)/w{k}(ij))/4; % why factor 4?
                Sr(rsk,csk) = sqw*Da{n,k}(lm,ij)*Mc(n,k);
                csk = csk + 1;                
            end            
        end
        rsk = rsk + 1;
    end
end


Ss = zeros(Np);
rsk = 1;
for n = 1:N+1
    for lm = 1:n^2
        csk = 1;
        for k = 1:N+1
            for ij = 1:k^2
                sqw = sqrt(w{n}(lm)/w{k}(ij))/4; % why factor 4?
                Ss(rsk,csk) = sqw*Db{n,k}(lm,ij)*Mc(n,k);
                csk = csk + 1;                
            end            
        end
        rsk = rsk + 1;
    end
end

St = zeros(Np);
rsk = 1;
for n = 1:N+1
    for lm = 1:n^2
        csk = 1;
        for k = 1:N+1
            for ij = 1:k^2
                sqw = sqrt(w{n}(lm)/w{k}(ij))/4; % why factor 4?
                Sra = .5*(1+a{n}(lm))*Sr(rsk,csk);
                Ssb = .5*(1+b{n}(lm))*Ss(rsk,csk);
                St(rsk,csk) = Sra + Ssb + sqw*VDMab{n,k}(lm,ij)*S1Dc(n,k);
                csk = csk + 1;                
            end            
        end
        rsk = rsk + 1;
    end
end

