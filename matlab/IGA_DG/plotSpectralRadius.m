NB = 3;
Ksub = 8;

re = linspace(-1,1,Ksub+1);
rp = linspace(-1,1,1000)';
Vp = bsplineVDM(NB,Ksub,rp);
plot(rp,Vp)
hold on
plot(re,0*re,'o')

%% p-ref

clear
Nvec = 1:7;
K1D = 1;

sk = 1;
for NB = Nvec;
    Ksub(sk) = NB;
%     rho(sk) = advec1D_spectra(NB,Ksub(sk),K1D,50);
    rho(sk) = Wave1D_spectra(NB,Ksub(sk),K1D,0);
%     rho(sk) = IPDG_IGA1D_spectra(NB,Ksub,K1D,0);
    sk = sk + 1;
end
plot(Nvec+Ksub,rho,'o--')
hold on

%% p-ref, polynomials

tau = 2;
Nvec = 1:8;
sk = 1;
for NB = Nvec;
    smoothKnots = 0;
    K1D = 1;
    Ksub = 1;
    rho(sk) = advec1D_spectra(NB,Ksub,K1D,smoothKnots,tau);
%     rho(sk) = Wave1D_spectra(NB,Ksub,K1D,smoothKnots,0);
%     rho(sk) = IPDG_IGA1D_spectra(NB,Ksub,K1D,0);
    sk = sk + 1;
end
plot(Nvec+1,rho,'x--')
hold on

%% href

clear
NB = 5;
Kvec = 32:32:256;
% Kvec = 2.^(2:7);
K1D = 1;
smoothKnots = 0;
tau = .5;
sk = 1;
for Ksub = Kvec
    Ksub
    [rho lam] = advec1D_spectra(NB,Ksub,K1D,smoothKnots,tau);
    rhoB(sk) = rho;
%     rhoB(sk) = Wave1D_spectra(NB,Ksub,K1D,smoothKnots,tau);
%     rhoB(sk) = IPDG_IGA1D_spectra(NB,Ksub,K1D,0);
    sk = sk + 1;
end
plot(Kvec,rhoB,'o--')
hold on
print_pgf_coordinates(Kvec,rhoB)

%%
% polynomials
tau = 0;
Ksub = 1;
sk = 1;
for K1D = Kvec
    rhoDG(sk) = advec1D_spectra(NB,Ksub,K1D,tau);
%     rhoDG(sk) = Wave1D_spectra(NB,Ksub,K1D,tau);
%     rho(sk) = IPDG_IGA1D_spectra(NB,Ksub,K1D,0);
    sk = sk + 1;
end
plot(Kvec,rhoDG,'x--')
hold on

% plot(Kvec,2*Kvec,'*--')
% hold on

% Nvec = [4:16];
% K1D = 1;
% Ksub = 1;
% clear rho
% sk = 1;
% for NB = Nvec
%     rho(sk) = advec1D_spectra(NB,Ksub,K1D,0);
%     %rho(sk) = Wave1D_spectra(NB,Ksub,K1D,0);
% %     rho(sk) = IPDG_IGA1D_spectra(NB,Ksub,K1D,0);
%     sk = sk + 1;
% end
% plot((Nvec+1),rho,'s--')
% hold on