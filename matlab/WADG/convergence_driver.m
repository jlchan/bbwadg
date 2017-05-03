clear

N = 4;
for ref = 1:4
    K1D = 2.^ref;    
    L2err_orig(ref,1) = Wave2D_RK(N,K1D,1);
    L2err_cmass(ref,1) = Wave2D_RK(N,K1D,0);
end

h = 2*.5.^(1:length(L2err_orig)); 
h  = h(:); 

%%

print_pgf_coordinates(h,L2err_orig)
print_pgf_coordinates(h,L2err_cmass)

loglog(h,L2err_orig,'o-')
hold on
loglog(h,L2err_cmass,'s-')

ids = 2:ref;
C1 = [h(ids).^0 log(h(ids))]\log(L2err_orig(ids));
C2 = [h(ids).^0 log(h(ids))]\log(L2err_cmass(ids));
title(sprintf('L2 convergence rate orig = %f, cmass = %f, expected rate %f',C1(2),C2(2),N+1))

format shorte
fprintf('%4.4e & ',L2err_orig);fprintf('\n')
fprintf('%4.4e & ',L2err_cmass);fprintf('\n')
format
rates = [C1(2), C2(2)]
