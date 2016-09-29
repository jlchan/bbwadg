clear
for N = 1:9
    
    [nnzs conds] = Advec2D_precon(N);
    nnzP(N)=nnzs(1);
    nnzPavg(N)=nnzs(2);
    condA(N) = conds(1);
    condDense(N) = conds(2);
    condSparse(N) = conds(3);
end

hold on
plot(condA,'o--')
plot(condDense,'s--')
plot(condSparse,'x--')