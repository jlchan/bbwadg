clear

for N = 1:4
    [r s t w] = hex_cubature(N);
    [V Dr Ds Dt] = hex_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CH(N) = max(abs(eig(K)));
    CHr(N) = max(abs(eig(Kr)));
    CHs(N) = max(abs(eig(Ks)));
    CHt(N) = max(abs(eig(Kt)));    
    
    [r s t w] = wedge_cubature(N);
    [V Dr Ds Dt] = wedge_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CW(N) = max(abs(eig(K)));
    CWr(N) = max(abs(eig(Kr)));
    CWs(N) = max(abs(eig(Ks)));
    CWt(N) = max(abs(eig(Kt)));   
    
    [r s t w] = pyr_cubature(N);
    [V Dr Ds Dt] = pyr_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CP(N) = max(abs(eig(K)));
    CPr(N) = max(abs(eig(Kr)));
    CPs(N) = max(abs(eig(Ks)));
    CPt(N) = max(abs(eig(Kt)));
    
    [r s t w] = tet_cubature(2*N);
    [V Dr Ds Dt] = tet_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CT(N) = max(abs(eig(K)));
    CTr(N) = max(abs(eig(Kr)));
    CTs(N) = max(abs(eig(Ks)));
    CTt(N) = max(abs(eig(Kt)));
    N
end

NN = 1:N;
%%
% pH = polyfit(log(NN),log(CHr),1)
% pW = polyfit(log(NN),log(CWr),1)
% pP = polyfit(log(NN),log(CPt),1)
% pT = polyfit(log(NN),log(CTr),1)

% mean([log(CHr(N)/CHr(N-1))/log(N/(N-1)) log(CHr(N-1)/CHr(N-2))/log((N-1)/(N-2))])
% mean([log(CWr(N)/CWr(N-1))/log(N/(N-1)) log(CWr(N-1)/CWr(N-2))/log((N-1)/(N-2))])
% mean([log(CPr(N)/CPr(N-1))/log(N/(N-1)) log(CPr(N-1)/CPr(N-2))/log((N-1)/(N-2))])
% mean([log(CTr(N)/CTr(N-1))/log(N/(N-1)) log(CTr(N-1)/CTr(N-2))/log((N-1)/(N-2))])

% plot(log(NN),log(CHr),'b.-','linewidth',2);hold on
% plot(log(NN),pH(1)*log(NN)+pH(2),'r.--')

% plot(log(NN),log(CWr),'b.-','linewidth',2);hold on
% plot(log(NN),pW(1)*log(NN)+pW(2),'r.--')

%%

if 0
    loglog(NN,CWr,'ro-','linewidth',2);hold on
    
    % loglog(NN,CPt,'ms-','linewidth',2);hold on
    loglog(NN,CTr,'k*-','linewidth',2);hold on
    
    hold on
    p = polyfit(NN,CHr,3);
    % plot(NN,p(1)*NN.^2 + p(2)*NN + p(3),'*--')
    % plot(NN,p(1)*NN.^3 + p(2)*NN.^2 + p(3)*NN + p(4),'*--')
    plot(NN,NN.^2,'bo--')
    plot(NN,NN.^3,'ro--')
    plot(NN,NN.^4,'ko--')
    % legend('r','s','t','Trace');set(gca,'fontsize',14)
    
    return
end
%%

% subplot(2,2,1)

semilogy(NN,sqrt(CH),'b.-','linewidth',2);hold on
plot(NN,2*(NN+1).*(NN+2)/2,'ro--','linewidth',2);hold on
legend('Sqrt Markov','2 x Trace');set(gca,'fontsize',14)
xlabel('N','fontsize',14)
% title('Hexahedra','fontsize',15)
% print(gcf,'-depsc','../docs/figs/markovTraceH.eps')
return
clf
% subplot(2,2,2)
% figure
semilogy(NN,sqrt(CW),'b.-','linewidth',2);hold on
plot(NN,2*(NN+1).*(NN+2),'ro--','linewidth',2);hold on
legend('Sqrt Markov','2 x Trace');set(gca,'fontsize',14)
xlabel('N','fontsize',14)
% title('Prism','fontsize',15)
% print(gcf,'-depsc','../docs/figs/markovTraceW.eps')

clf
% subplot(2,2,3)
% figure
semilogy(NN,sqrt(CP),'b.-','linewidth',2);hold on
plot(NN,2*(NN+1).*(NN+3),'ro--','linewidth',2);hold on
legend('Sqrt Markov','2 x Trace');set(gca,'fontsize',14)
xlabel('N','fontsize',14)
% title('Pyramid','fontsize',15)
% print(gcf,'-depsc','../docs/figs/markovTraceP.eps')

clf
% subplot(2,2,4)
% figure
semilogy(NN,sqrt(CT),'b.-','linewidth',2);hold on
plot(NN,2*(NN+1).*(NN+3)/2,'ro--','linewidth',2);hold on
legend('Sqrt Markov','2 x Trace');set(gca,'fontsize',14)
xlabel('N','fontsize',14)
% title('Tetrahedra','fontsize',15)
% print(gcf,'-depsc','../docs/figs/markovTraceT.eps')
return
