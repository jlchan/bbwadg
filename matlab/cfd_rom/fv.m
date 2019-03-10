clear
K = 100;
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

for e = 2:K-1 
    Sk{e} = zeros(size(S));
    Sk{e}(e-1:e+1,e-1:e+1) = .5*[0 1 0;-1 0 1;0 -1 0];   
end
Sk{1} = zeros(size(S)); 
Sk{1}(1:2,1:2) = .5*[0 1;-1 0];
Sk{1}(1,K) = -1/2; Sk{1}(K,1) = 1/2;
Sk{K}(K-1:K,K-1:K) = .5*[0 1;-1 0];
Sk{K}(1,K) = -1/2; Sk{K}(K,1) = 1/2;

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
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);
% fD = @(uL,uR) 0*max(abs(uL),abs(uR)).*(uR-uL);

% fS = @(uL,uR) (uL + uR)/2;
% fD = @(uL,uR) (uR-uL);

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
     
u = 1-sin(pi*x);

% make reduce basis
if 1
    N = ceil(K/2);
    %Vrfull = orth([x(:).^(0:N) randn(K,2)]);
    Vrfull = orth([x(:).^(0:N)]);
    
    Nmodes = ceil(size(Vrfull,2)/2);    
    Vr = Vrfull(:,1:Nmodes);
    Vrp = Vr;
else
    
    load Usnap_burgers
    [Ur,Sr,Vr] = svd(Usnap,0);
    
    Nmodes = 5;
    Vrfull = [x(:).^0 Ur(:,1:2*Nmodes)];
    Vr = orth([x(:).^0 Vrfull(:,1:Nmodes)]);
    Vrp = Vr;
    
end

opt = 0;
if opt==0
    
%     plot(x,Vr*Vr'*u,'bo')
%     hold on;
%     plot(x,u,'r--')
%     return

    % project
    u = Vr'*u;    
    
    SD = abs(S);
%     S = Vr*Vr'*S*Vr*Vr';
    invM = eye(Nmodes)/dx;
    
elseif opt==1
    
    % project
    u = Vr'*u;
    
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
            
elseif opt==2        
    Vr = eye(K);
    Vrp = Vr;
    SD = abs(S);
    invM = eye(K)/dx;
end


dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));

for i = 1:Nsteps
    for INTRK = 1:5                        
        [ux uy] = meshgrid(Vr*u);
        rhs = -Vr'*(sum(S.*fS(ux,uy),2) + sum(abs(SD).*fD(ux,uy),2));        
%         rhs = -Vr'*(sum(S.*fS(ux,uy),2));  
        
        if (INTRK==5)
            rhstest(i) = u'*rhs;
        end
        rhs = invM*rhs;

        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res;                
    end        
    
    if mod(i,1) == 0
        plot(x,Vrp*u,'o--')
        axis([-1,1,-2,2])
        title(sprintf('max rhstest = %g\n',max(abs(rhstest(i)))))
        drawnow
    end
end
