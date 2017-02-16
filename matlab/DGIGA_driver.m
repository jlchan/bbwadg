%% 1D
clear
sk = 1;
for NB = 2:9
    K1D = 9; Ksub = 1; % DG
    K1D = 4; Ksub = NB; % DG-IGA
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
    sk = sk + 1;
    NB
end
'done'

loglog(work,L2err,'x--')
hold on

%% 2D

clear
sk = 1;
for NB = 2:9
    K1D = 7; Ksub = 1; % DG
    K1D = 3; Ksub = NB; % DG-IGA
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
    NB
end
'done'

loglog(dofs,L2err,'x--','linewidth',2)
xlabel('Total work','fontsize',14)
ylabel('L2 error','fontsize',14)
hold on