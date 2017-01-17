clear
hybridgGlobals3D
hybridgGlobalFlags 

useSEM = 1;
for N = 1:9
    
    %% hex
    Nfaces= 6; K = 1;
    
    VX = [-1  1  1 -1 -1  1  1 -1]'; VY = [-1 -1  1  1 -1 -1  1  1]'; VZ = [-1 -1 -1 -1  1  1  1  1]';
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    hybridgStartUp
            
    Vf = hex_basis(N,rfH,sfH,tfH); % eval map @ cubature
    Mf = Vf'*spdiag(wsJ(1:NfcH,1))*Vf; % incorporate face scalings
    V = hex_basis(N,rqH,sqH,tqH);    
    M = V'*spdiag(wJ(1:NcH,1))*V;        
    CH(N) = max(abs(eig(M\Mf)));
     
    [r s t w] = hex_cubature(N);
    [V Dr Ds Dt] = hex_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CMH(N) = max(abs(eig(M\K)));   
    
    
    %% wedge
    Nfaces = 5; K = 1;
    u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
    VX = v(:); VY = w(:); VZ = u(:); % flipping coordinates for Gmsh
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
        
    hybridgStartUp
    
    V = wedge_basis(N,rfW,sfW,tfW); % eval map @ cubature
    Mf = V'*spdiag(wsJ(1:NfcW,1))*V;
    V = wedge_basis(N,rqW,sqW,tqW);    
    M = V'*spdiag(wJ(1:NcW,1))*V;    
    CW(N) = max(abs(eig(M\Mf)));
    
    [r s t w] = wedge_cubature(N);
    [V Dr Ds Dt] = wedge_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CMW(N) = max(abs(eig(M\K)));    
    
    
    %% pyramid
    Nfaces = 5; K = 1;
    VX = [ -1   1   1  -1  -1 ]';VY = [ -1  -1   1   1  -1 ]'; VZ = [ -1  -1  -1  -1   1 ]';
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    hybridgStartUp
    
    V = pyr_basis(N,rfP,sfP,tfP); % eval map @ cubature
    Mf = V'*spdiag(wsJ(1:NfcP,1))*V;
    V = pyr_basis(N,rqP,sqP,tqP);    
    M = V'*spdiag(wJ(1:NcP,1))*V;    
    
    CP(N) = max(abs(eig(M\Mf)));
    
    [r s t w] = pyr_cubature(N);
    [V Dr Ds Dt] = pyr_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CMP(N) = max(abs(eig(M\K)));    
        
    %% tet
    Nfaces = 4; K = 1;
    VX = [-1  1 -1 -1 ]'; VY = [-1 -1  1 -1 ]'; VZ = [-1 -1 -1  1 ]';
    
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    hybridgStartUp
    
    V = tet_basis(N,rfT,sfT,tfT); % eval map @ cubature
    Mf = V'*spdiag(wsJ(1:NfcT,1))*V;
    V = tet_basis(N,rqT,sqT,tqT);    
    M = V'*spdiag(wJ(1:NcT,1))*V;    
    
    CT(N) = max(abs(eig(Mf)));
    
    [r s t w] = tet_cubature(2*N);
    [V Dr Ds Dt] = tet_basis(N,r,s,t); % eval map @ cubature    
    Kr = Dr'*diag(w)*Dr; Ks = Ds'*diag(w)*Ds; Kt = Dt'*diag(w)*Dt;
    K = Kr + Ks + Kt;
    CMT(N) = max(abs(eig(M\K))); 
    
    %%
    disp(sprintf('done with N = %d',N))
end

NN = 1:N;
% CH
% CMH

keyboard
%% plot hex
figure
plot(NN,.5*CH,'b.-','linewidth',2) 
hold on
plot(NN,sqrt(CMH),'ro--','linewidth',2)
set(gca,'fontsize',14)
legend('Trace constant','Sqrt Markov constant')
% print(gcf,'-depsc','../docs/figs/hexTrace.eps')
%% wedge
figure
plot(NN,.5*CW,'b.-','linewidth',2)
hold on
plot(NN,sqrt(CMW) ,'ro--','linewidth',2)
set(gca,'fontsize',14)
legend('Trace constant','Sqrt Markov constant')
% print(gcf,'-depsc','../docs/figs/wedgeTrace.eps')
%% pyramid
figure
plot(NN,.5*CP,'b.-','linewidth',2)
hold on
plot(NN,sqrt(CMP),'ro--','linewidth',2)
set(gca,'fontsize',14)
legend('Trace constant','Sqrt Markov constant')
% print(gcf,'-depsc','../docs/figs/pyrTrace.eps')
%% tet
figure
plot(NN,.5*CT,'b.-','linewidth',2)
hold on
plot(NN,sqrt(CMT),'ro--','linewidth',2)
legend('Trace constant','Sqrt Markov constant')
set(gca,'fontsize',14)
% print(gcf,'-depsc','../docs/figs/tetTrace.eps')

