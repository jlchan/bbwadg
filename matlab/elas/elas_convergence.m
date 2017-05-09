clc
% clear

% L2 errors
Kvec = [2 4 8 16];
for N = 5;1:5
    for i = 1:length(Kvec);
        K1D = Kvec(i);
        h(i) = 1/K1D;
        disp(sprintf('======== on N = %d, K1D = %d\n',N,K1D));
        
        harmonic_err{N}(i) = Elasticity2D(N,K1D) % harmonic sol
        
        rayleigh_err{N}(i) = Elasticity2D_RayleighWave(N,K1D) % mu = 1
        
        lamb_err{N}(i) = Elasticity2D_Lamb(N,K1D)        

        stonely_err{N}(i) = Elasticity2D_Stonely(N,K1D)
        
%         ref_err{N}(i) = ElasticityReferenceSolution2D(N,K1D);                
    end            
end
%%
for N = 1:5
%     err = harmonic_err{N};
%     err = rayleigh_err{N};    
%     err = lamb_err{N};
    err = stonely_err{N};
    
    loglog(h,err,'o--')
    hold on
    r = N+1/2;
    scale = err(1)/h(1)^(r);
    loglog(h,scale*h.^(r),'k--')
    C(N) = (log(err(end))-log(err(end-1))) / (log(h(end))-log(h(end-1)));
    a = 1;
    C2=[ones(a+1,1) log(h(end-a:end))']\log(err(end-a:end)'); C(N) = C2(2);
    print_pgf_coordinates(h,err)    
end
C

%%
% clc

for N=1:6
%     err = harmonic_err{N};
    %err = rayleigh_err{N};
%     err = lamb_err{N};
%     err = stonely_err{N};
    err = ref_err{N};
    print_pgf_coordinates(h,err)    
    a = 1;
    C = [ones(a+1,1) log(h(end-a:end))']\log(err(end-a:end)'); disp(C(2))
end
%%

Kvec2 = 4;
% Kvec2 = [4 8];
muvec = [1 .1 .01 .001 .0001]; lambdavec = ones(size(muvec));
% lambdavec = [1 2 4 8 16]; muvec = ones(size(muvec)); % other approach

for N = 2:4
    for i = 1:length(Kvec2);
        for j = 1:length(muvec)
            mu = muvec(j);
            lambda = lambdavec(j);
            K1D = Kvec2(i);
            disp(sprintf('======== on N = %d, K1D = %d, mu = %f, lam = %f\n',N,K1D,mu,lambda));
            [l2err errU errS] = Elasticity2D_RayleighWave(N,K1D,mu,lambda); % mu = 1
            rayleigh_incomp_errU{i}(N,j) = errU;
            rayleigh_incomp_errS{i}(N,j) = errS;            
        end
    end
end

% herr = [0.5,0.1732;0.25,0.0439;0.125,0.0110;0.0625,0.0027];
% h = herr(:,1); err = herr(:,2);
% C = [ones(size(h)) log(h)]\log(err(:)); C(2)
