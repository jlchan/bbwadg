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
fD = @(uL,uR) (uR-uL);
% fD = @(uL,uR) 0*max(abs(uL),abs(uR)).*(uR-uL);

% fS = @(uL,uR) (uL + uR)/2;
% fD = @(uL,uR) (uR-uL);

SD = abs(S);

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

interval = 1;
dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));

Usnap = zeros(K,ceil(Nsteps/interval));
sk = 1;
for i = 1:Nsteps
    for INTRK = 1:5                        
        [ux uy] = meshgrid(u);
%         rhs = -(sum(S.*fS(ux,uy),2) + sum(abs(SD).*fD(ux,uy),2));        
        rhs = -(sum(S.*fS(ux,uy),2) + max(abs(u(:)))*sum(abs(SD).*fD(ux,uy),2));        
%         rhs = -Vr'*(sum(S.*fS(ux,uy),2));                  
        rhs = rhs/dx;

        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res;                
    end        
    if mod(i,interval)==0
        Usnap(:,sk) = u;
        sk = sk + 1;
    end
    
    if mod(i,5) == 0
        plot(x,u,'o--')
        axis([-1,1,-2,2])       
        drawnow
    end
end
