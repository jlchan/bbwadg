
Kvec = [8 16 32 64];
L2errGQ = {};
L2errGLL = {};
for N = 1:4    
    for K1D = Kvec
        fprintf('N = %d, K1D = %d\n',N,K1D);
        
        CFL = .4;
        [rq1D wq1D] = JacobiGQ(0,0,N);
        err = Euler2D_skew(N,K1D,rq1D,wq1D,CFL);
        L2errGQ{N,K1D} = err;
        
%         CFL = .25;
%         [rq1D wq1D] = JacobiGL(0,0,N);
%         err = Euler2D_skew(N,K1D,rq1D,wq1D,CFL);
%         L2errGLL{N,K1D} = err;
        
    end
end