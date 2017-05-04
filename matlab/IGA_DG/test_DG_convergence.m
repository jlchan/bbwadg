clc
clear
smoothKnots = 75;
for NB = 2:5;    
    e0 = []; 
    ndofs = [];
    for K1D = [2 4 8]
        Kvec = 2.^(2:7);
        Kvec = Kvec(Kvec*K1D < 2^7);
        clear L2err dofs
        sk = 1;
        for Ksub = Kvec
%             [L2e nndofs] = IPDG_IGA1D(NB,Ksub,K1D,smoothKnots);
            [L2e nndofs] = Wave1D_IGA(NB,Ksub,K1D,smoothKnots);
            L2err(sk) = L2e;
            dofs(sk) = nndofs;
            sk = sk + 1;
        end
        if NB==2
            semilogy(dofs,L2err,'o--','linewidth',2)
        elseif NB==3
            loglog(dofs,L2err,'x--','linewidth',2)
        elseif NB==4
            loglog(dofs,L2err,'^--','linewidth',2)
        elseif NB==5
            loglog(dofs,L2err,'*--','linewidth',2)
        end
        
        hold on
        
        e0 = [e0 L2err(1)];
        ndofs = [ndofs dofs(1)];

        print_pgf_coordinates(dofs,L2err)
        
    end
%     C = polyfit(log(ndofs(2:end)),log(e0(2:end)),1);
%     loglog(ndofs,exp(C(1)*log(ndofs) + C(2)),'k-','linewidth',2)
%     C(1)
%     disp('% === ')
%     print_pgf_coordinates(ndofs,exp(C(1)*log(ndofs) + C(2)))
%     disp('% === ')
    
end
grid on 
set(gca,'fontsize',14)
xlabel('Number of degrees of freedom','fontsize',14)

%% 2D

clc
clear
smoothKnots = 0;
for NB = 2:5;    
    e0 = []; 
    ndofs = [];
    for K1D = [2 4 8]
        Kvec = 2.^(1:4);
        Kvec = Kvec(Kvec*K1D < 2^7);
        clear L2err dofs
        sk = 1;
        for Ksub = Kvec
            [L2e nndofs] = Wave2D_IGA(NB,Ksub,K1D,smoothKnots);
%             [L2e nndofs] = IPDG_IGA1D(NB,Ksub,K1D,smoothKnots);
%             [L2e nndofs] = Wave1D_IGA(NB,Ksub,K1D,smoothKnots);
            L2err(sk) = L2e;
            dofs(sk) = nndofs;
            sk = sk + 1;
        end
        if NB==2
            semilogy(dofs,L2err,'o--','linewidth',2)
        elseif NB==3
            loglog(dofs,L2err,'x--','linewidth',2)
        elseif NB==4
            loglog(dofs,L2err,'^--','linewidth',2)
        elseif NB==5
            loglog(dofs,L2err,'*--','linewidth',2)
        end
        
        hold on
               
        print_pgf_coordinates(dofs,L2err)
        
    end
    
end
grid on 
set(gca,'fontsize',14)
xlabel('Number of degrees of freedom','fontsize',14)