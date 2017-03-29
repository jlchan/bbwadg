%% 1D
clear
sk = 1;
for NB = 10:20
    K1D = 9; Ksub = 1; % DG
    K1D = 4; Ksub = 1; % p-ref
%     K1D = 4; Ksub = NB; % DG-IGA
    
    K1D = 8; Ksub = 1; % DG-IGA
%     K1D = 8; Ksub = NB; % DG-IGA
    if Ksub==1
        dt = .5 / (NB^2*K1D);
    else
        dt = .5 / (NB*Ksub*K1D);
    end

    L2err(sk) = advec1D_IGA(NB,Ksub,K1D,dt);
    N = NB+Ksub-1;
    Nsteps = ceil(1/dt);
    d = 1;
    work(sk) = Nsteps*(N+1)^(d+1)*K1D^2;  % for 3D
    dofs(sk) = (N+1)*K1D;
    sk = sk + 1;
    NB
end
'done'

loglog(dofs,L2err,'x--')
hold on

%% more 1D
clear

NB = 6;
for Ksub = 1:8;
    sk = 1;
    
    K1D = 1;
    for i = 1:4
        %     K1D = 9; Ksub = 1; % DG
        K1D = 2*K1D; %Ksub = NB;
        dt = .5 / (NB^2*K1D);
        L2err(sk) = advec1D_IGA(NB,Ksub,K1D,dt);
        N = NB+Ksub-1;
        dofs(sk) = (N+1)*K1D;
        sk = sk + 1;
    end
    %loglog(dofs,L2err,'x--','markersize',8,'linewidth',2,'DisplayName',sprintf('DG, p = %d',NB))
    loglog(dofs,L2err,'x--','markersize',8,'linewidth',2,'DisplayName',sprintf('p = %d, Ksub = %d',NB,Ksub))
    hold on
    NB
end
legend('show','Location','Best')
set(gca,'fontsize',15)
grid on
title('h-refinement','fontsize',15)
axis([4 500  1e-8  10])


%% p-ref in 1D

clear
for Ksubb = [1 2 4 8]
    sk = 1;    
    for NB = 3:9

        if Ksubb==1
            K1D = 16;
        elseif Ksubb==2
            K1D = 15;
        elseif Ksubb==4
            K1D = 14;
        else
            K1D = 13;
        end
        Ksub = 1;
        dt = .5 / (NB^2*K1D);
        L2err(sk) = advec1D_IGA(NB,Ksub,K1D,dt);
        N = NB+Ksub-1;
        dofs(sk) = (N+1)*K1D;
        sk = sk + 1;
    end
    
    loglog(dofs,L2err,'x--','markersize',8,'linewidth',2,'DisplayName',sprintf('Ksub = %d, K1D = %d',Ksub,K1D))
    hold on
end

set(gca,'fontsize',15)
grid on
title('p-refinement (p = 3 to p = 7)','fontsize',15)
axis([4 500  1e-8  10])

legend('show','Location','Best')


%% 2D

clear
sk = 1;
K1D = 1/2;
for NB = 2:5
%     K1D = 7; Ksub = 1; % DG
%     K1D = 3; Ksub = NB; % DG-IGA
%     K1D = 2; Ksub = 1; % DG
    K1D = 2; Ksub = 2*NB; % DG-IGA
    disp(sprintf('NB = %d\n',NB))
    
    if Ksub==1
        dt = .5 / (NB^2*K1D);
    else
        dt = .5 / (NB*Ksub*K1D);
    end

    L2err(sk) = WaveQuad_IGA(NB,Ksub,K1D,dt);
    
    N = NB+Ksub-1;
    
    Nsteps = ceil(1/dt);
    d = 2;
    if Ksub==1
        work(sk) = Nsteps*(N+1)^(d+1)*K1D^(2*d);  % for SEM-DG - one square matvec
    else
        MultCost = 4*(N+1)^d * (2*N-1); % estimate quadrature needs 2*N-1 nodes?
        
        work(sk) = Nsteps*MultCost*K1D^(2*d);  % for IGA-DG - factor in more matmults (3 quadrature matmults + larger dim?)
    end
    dofs(sk) = (N+1)^(d)*K1D^(d);
    sk = sk + 1;    
end
'done'

loglog(dofs,abs(L2err),'x--','linewidth',2)
%xlabel('Total work','fontsize',14)
xlabel('Dofs','fontsize',14)
ylabel('L2 error','fontsize',14)
hold on