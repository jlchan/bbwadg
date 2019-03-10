clear
K = 200;
FinalTime = 2;
xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 2/K;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

SD = abs(S);
Kr = 2*eye(size(S))-abs(S);

% w = ones(K,1);
% w = [1; round(rand(K-2,1)); 1];
% id = find(w);
% Sr = 0;    
% C = zeros(K);
% for e = 1:K
%     Sr = Sr + w(e)*Sk{e};
%     C(:,e) = sum(Sk{e},2);
% end

fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
% fS = @(uL,uR) (uL.^2 + uR.^2)/4;
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);
% fD = @(uL,uR) (uR-uL);
% fD = @(uL,uR) 0*max(abs(uL),abs(uR)).*(uR-uL);

% fS = @(uL,uR) (uL + uR)/2;
% fD = @(uL,uR) 0*(uR-uL);

%% 

rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];
     
u = -sin(pi*x);

opt = 2;

% make reduce basis
if opt==0 % FOM        
    Nmodes = K;
    Vr = eye(K);
    Vrp = Vr;    
    invM = eye(K)/dx;    
elseif opt==1
    
    N = ceil(K/2);
    %Vrfull = orth([x(:).^(0:N) randn(K,2)]);
    Vrfull = orth([x(:).^(1:N)]);
%     Vrfull = [x(:).^0 orth([x(:).^(1:N)])];
%     Vrfull = [x(:).^(0:4)];
    
    Nmodes = ceil(size(Vrfull,2)/2);    
    Vr = Vrfull(:,1:Nmodes);
    Vrp = Vr;
       
    invM = inv(dx*Vr'*Vr);
    Pr = invM*Vr'*dx;
    
    % project solution
    u = Pr*u; 
    
else    
%     load Usnap_burgers_moving
    load Usnap_burgers_GLF
%     load Usnap_burgers
    [Ur,Sr,Vr] = svd(Usnap,0);
    
    Nmodes = 35;
    Vr = [x(:).^0 Ur(:,1:Nmodes)];    
    
    Vrp = Vr;
    
    invM = inv(dx*Vr'*Vr);
    Pr = invM*Vr'*dx;
    
    % project solution
    u = Pr*u;            
    
end



opt = 0;
if opt==0
    
    % no hyperreduction
    
    % viscosity matrix
    Kr = Vr'*Kr*Vr;
    
elseif opt==1        
    
    Sfull = S;
    [Utest,Stest,~] = svd([Vr S*Vr],0); % enrich with "least squares" basis 
    Vtest = Utest(:,diag(Stest) > 1e-12); % include nonzero energy modes
    
%     % compress 
%     Prr = (Vtest'*Vtest)\(Vtest');    
%     Stest = Prr'*(Vtest'*Sfull*Vtest)*Prr;
%     SDtest = abs(Sfull);
    
    % viscosity matrix
    Kr = Vr'*Kr*Vr;
    
    id = get_DEIM_ids(Vtest);    
    Vrred = Vtest(id,:);    
    Stest = Vrred'\((Vtest'*Sfull*Vtest)/Vrred); 
    SDtest = abs(Sfull(id,id)); 
    Vr = Vr(id,:);
    
    % hyperreduction
    S = Stest; 
    SD = SDtest; 
    
    
elseif opt==2
    
    % hyperreduction with deim
    id = get_DEIM_ids(Vrfull);
    Vrred = Vr(id,:);
    w = Vrfull(id,:)'\(sum(dx*Vrfull',2)); % quad weights with constraint on const exactness
    
    %     WL = pinv([2*Vrfull(id,:)*Vrfull(id,:)' ones(size(id))';
    %         ones(size(id)) 0])*[2*sum(dx*Vrfull',2); dx*K];
    %     w = WL(1:end-1);
    
    plot(x(id),id*0,'o')
    
    return
    
    M = (Vrred'*diag(w)*Vrred);
    Pr = M\(Vrred'*diag(w));
    Sfull = S;
    S = Pr'*(Vr'*Sfull*Vr)*Pr;
    SD = abs(Sfull(id,id));
    %     SD = abs(S);
    
    Vr = Vrred;
        
    invM = inv(Vrred'*diag(w)*Vrred);
           
end

dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));

for i = 1:Nsteps
    for INTRK = 1:5          
        uq = Vr*u;
        [ux uy] = meshgrid(uq);        
%         rhs = Vr'*(sum(S.*fS(ux,uy),2) + sum(abs(SD).*fD(ux,uy),2)); 
        rhs = Vr'*(sum(S.*fS(ux,uy),2)); 
%         rhs = Vr'*(sum(S.*fS(ux,uy),2));                            

        if (INTRK==5)
            rhstest(i) = u'*rhs;
        end
        
        tau = 1;
        rhs = -invM*(rhs + tau*max(abs(uq(:)))*Kr*u); % invert mass + add art visc

        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res;                
    end        
    
    if mod(i,5) == 0     
%         plot(diag(S.*fS(ux,uy),1))
%         plot(diag(abs(SD).*fD(ux,uy),1),'o')
        plot(x,Vrp*u,'bo--')
        hold on
        plot(x,Usnap(:,i),'rx-')
        hold off
        axis([-1,1,-1,2.5])
        title(sprintf('max rhstest = %g\n',max(abs(rhstest(i)))))
        drawnow
    end
end
