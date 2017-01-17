clear
hybridgGlobals3D
hybridgGlobalFlags
useSEM=1; useLSC=0;
a = .0;
    
for N = 1:9
    
    %% hex
    Nfaces= 6; K = 1;
    
    VX = [-1  1  1 -1 -1  1  1 -1]'; VY = [-1 -1  1  1 -1 -1  1  1]'; VZ = [-1 -1 -1 -1  1  1  1  1]';
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    VX = VX + a*randn(size(VX));
    VY = VY + a*randn(size(VY));
    VZ = VZ + a*randn(size(VZ));
    
%     VX = VX*2;  VY = VY*2; VZ = VZ*2;
    hybridgStartUp
            
    Vf = hex_basis(N,rfH,sfH,tfH); % eval map @ cubature
    Mf = Vf'*spdiag(wsJ(1:NfcH,1))*Vf; % incorporate face scalings
    V = hex_basis(N,rqH,sqH,tqH);    
    M = V'*spdiag(wJ(1:NcH,1))*V;    
    CH(N) = max(abs(eig(M\Mf)));
        
%     dP = sum(wfH.*sJ(1:NfcH,1)); %24; 
%     P = sum(wJ(1:NcH,1)); %8;    
    dP = max(JsB(1:NfcH,1))*24; % Js = change in area from reference
    P = min(J(1:NcH,1))*8;     % J  = change in volume from reference
    
    bCH(N) = 4*(N+1).^2;%(N+1).*(N+1) * dP/P;        
        
    %% wedge
    Nfaces = 5; K = 1;
    u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
    VX = v(:); VY = w(:); VZ = u(:); % flipping coordinates for Gmsh
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    VX = VX + a*randn(size(VX));
    VY = VY + a*randn(size(VY));
    VZ = VZ + a*randn(size(VZ));
    hybridgStartUp
    
    Vf = wedge_basis(N,rfW,sfW,tfW); % eval map @ cubature
    Mf = Vf'*spdiag(wsJ(1:NfcW,1))*Vf;
    V = wedge_basis(N,rqW,sqW,tqW);    
    M = V'*spdiag(wJ(1:NcW,1))*V;    
%     keyboard
    CW(N) = max(abs(eig(M\Mf)));
    
%     dP = sum(wfW.*sJ(1:NfcW,1)); %4*sqrt(2)+12; 
%     P = sum(wJ(1:NcW,1)); %4;
    dP = max(JsB(1:NfcW,1))*(4*sqrt(2)+12); % Js = change in area from reference
    P = min(J(1:NcW,1))*4;                 % J  = change in volume from reference    
    bCW(N) = (sqrt(2)+3)*(N+1).*(N+2);

    %% pyramid
    Nfaces = 5; K = 1;
    VX = [ -1   1   1  -1  -1 ]';VY = [ -1  -1   1   1  -1 ]'; VZ = [ -1  -1  -1  -1   1 ]';
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    
    VX = VX + a*randn(size(VX)); 
    VY = VY + a*randn(size(VY)); 
    VZ = VZ + a*randn(size(VZ)); 
    hybridgStartUp
    
    Vf = pyr_basis(N,rfP,sfP,tfP); % eval map @ cubature
    Mf = Vf'*spdiag(wsJ(1:NfcP,1))*Vf;
    V = pyr_basis(N,rqP,sqP,tqP);    
    M = V'*spdiag(wJ(1:NcP,1))*V;    
    
    CP(N) = max(abs(eig(M\Mf)));
    
%     dP = sum(wfP.*sJ(1:NfcP,1)); %8+4*sqrt(2); 
%     P = sum(wJ(1:NcP,1)); %8/3;
    
    dP = max(JsB(1:NfcP,1))*(4*sqrt(2)+8); % Js = change in area from reference
    P = min(J(1:NcP,1))*(8/3);                 % J  = change in volume from reference    
    bCP(N) = (sqrt(2)+2)*(N+1).*(N+3);%(N+1).*(N+3) * dP/P;

    %% tet
    Nfaces = 4; K = 1;
    VX = [-1  1 -1 -1 ]'; VY = [-1 -1  1 -1 ]'; VZ = [-1 -1 -1  1 ]';
    
    EToV = 1:length(VX); EToE = zeros(1,Nfaces); EToF = zeros(1,Nfaces);
    VX = VX + a*randn(size(VX));
    VY = VY + a*randn(size(VY));
    VZ = VZ + a*randn(size(VZ));
        
    hybridgStartUp
    
    V = tet_basis(N,rfT,sfT,tfT); % eval map @ cubature
    Mf = V'*spdiag(wsJ(1:NfcT,1))*V;
    V = tet_basis(N,rqT,sqT,tqT);    
    M = V'*spdiag(wJ(1:NcT,1))*V;    
    
    CT(N) = max(abs(eig(Mf)));
       
    %%
    
%     dP = sum(wfT.*sJ(1:NfcT,1)); %2*sqrt(3) + 6; 
%     P = sum(wJ(1:NcT,1)); %4/3;    
    
    dP = max(JsB(1:NfcT,1))*(2*sqrt(3)+6);      % Js = change in area from reference
    P = min(J(1:NcT,1))*(4/3);                 % J  = change in volume from reference    
    bCT(N) = (sqrt(3)+3)/2 * (N+1)*(N+3);%(N+1).*(N+3)/2 * dP/P;
    
    %%
    disp(sprintf('done with N = %d',N))
end

NN = 1:N;

% plot(NN,bCH./CH,'b.-')
% hold on
% plot(NN,bCW./CW,'ro-')
% plot(NN,bCP./CP,'ms-')
% plot(NN,bCT./CT,'k^-')
keyboard
% return

%% plot hex
plot(NN,CH,'bs-','linewidth',2) 
hold on
plot(NN,bCH,'ro--','linewidth',2)
xlabel('N','fontsize',14)
axis([1 7 0 350])
set(gca,'fontsize',14)
legend('Computed constant','Bound')
print(gcf,'-depsc','../docs/figs/hexTrace.eps')
% wedge
clf
plot(NN,CW,'bs-','linewidth',2)
hold on
plot(NN,bCW ,'ro--','linewidth',2)
axis([1 7 0 350])
xlabel('N','fontsize',14)
set(gca,'fontsize',14)
legend('Computed constant','Bound')
print(gcf,'-depsc','../docs/figs/wedgeTrace.eps')
% pyramid
clf
plot(NN,CP,'bs-','linewidth',2)
hold on
plot(NN,bCP,'ro--','linewidth',2)
axis([1 7 0 350])
xlabel('N','fontsize',14)
set(gca,'fontsize',14)
legend('Computed constant','Bound')
print(gcf,'-depsc','../docs/figs/pyrTrace.eps')
% tet
clf
plot(NN,CT,'bs-','linewidth',2)
hold on
plot(NN,bCT,'ro--','linewidth',2)
axis([1 7 0 350])
xlabel('N','fontsize',14)
legend('Computed constant','Bound')
set(gca,'fontsize',14)
print(gcf,'-depsc','../docs/figs/tetTrace.eps')

keyboard
