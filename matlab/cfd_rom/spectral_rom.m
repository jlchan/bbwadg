% SBP spectral method

clear
N = 25;
x = linspace(-1,1,2*N+2)'; % interpolatory
x = x(1:end-1);
dx = x(2)-x(1);

sk = 1;
D = zeros(2*N+1);
for k = -N:N
    lamk = 1i*k*pi;
    Vq(:,sk) = exp(lamk*x)/sqrt(2);        
    D(sk,sk) = lamk;
    sk = sk + 1;
end
% plot(x,Vq)

W = dx*eye(2*N+1); % "weight" matrix
M = real(Vq'*W*Vq);
Pq = M\(Vq'*W);
D = real(Vq*D*Pq); % imag part = 0

K = (D'*W*D); %/dx;
Q = W*D; 

%% make 2D

Np = (2*N+1)^2;
[x y] = meshgrid(x);

%% flux

% Burgers
fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);

%% build reduced model

load Usnap_spectral

U0 = Usnap(:,1);
Us = Usnap;
[Vs, Ss,~] = svd(Us,0);
sig = diag(Ss);

Nmodes = 25;
tol = sqrt(sum(sig(Nmodes+1:end).^2)./sum(sig.^2));

Vr = Vs(:,1:Nmodes);
Vrp = Vr;

% full ops
Qxfull = kron(Q,W);
Qyfull = kron(W,Q);
Kfull = kron(K,W) + kron(W,K); % d2u/dx2 + d2u/dy2

% test space
[Vtest, Stest, ~] = svd([ones(Np,1) Vr Qxfull*Vr Qyfull*Vr],0);
sigt = diag(Stest);
sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
Vtest = Vtest(:,sigerr > tol);

u = Vr'*U0(:,1);
fprintf('initial err = %g\n',sqrt(sum(dx^2*(Vrp*u-U0).^2)))

% % % full ops
% Qx = Vtest*(Vtest'*Qxfull*Vtest)*Vtest';
% Qy = Vtest*(Vtest'*Qyfull*Vtest)*Vtest';
% K = Vtest*(Vtest'*Kfull*Vtest)*Vtest';
% invMVrT = (1/dx^2)*Vr'; % orth basis
% Pr = Vr';

% hyperreduc
if 1
    Vtarget = Vtest;
%     Vtarget = Vr;
    Vmass = zeros(Np,size(Vtarget,2)*(size(Vtarget,2)+1)/2);
    sk = 1;
    for i = 1:size(Vtarget,2)
        for j = i:size(Vtarget,2)
            Vmass(:,sk) = Vtarget(:,i).*Vtarget(:,j);
            sk = sk + 1;
        end
    end
    b = sum(Vmass'*dx^2,2);
    
    maxpts = Nmodes*10;
    id = get_empirical_cubature(Vmass,b,tol,maxpts);
    wr = Vmass(id,:)'\b;    
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));    
    Pr = Mr\(Vr(id,:)'*diag(wr));
    invMVrT = Mr\(Vr(id,:)');
    Vr = Vr(id,:);
    
    % make new ops
    Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
    Ptest = Mtest\(Vtest(id,:)'*diag(wr));    
    Qx = Ptest'*(Vtest'*Qxfull*Vtest)*Ptest;
    Qy = Ptest'*(Vtest'*Qyfull*Vtest)*Ptest;
    K = Ptest'*(Vtest'*Kfull*Vtest)*Ptest;
end


%%

FinalTime = 2;

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


Qxy = (Qx+Qy);

dt = .5*dx;
% dt = dt/2.5;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));
rhs = zeros(size(u));

for i = 1:Nsteps
    for INTRK = 1:5
        
        uq = Vr*u;
        [ux uy] = meshgrid(uq);
        rhs = sum(Qxy.*fS(ux,uy),2);
        
        if INTRK==5
            rhstest(i) = sum(sum(uq.*rhs));
        end
        
        tau = .1*dx;
        rhs = -invMVrT*(rhs + tau*K*uq);
        
%         if INTRK==5
%             rhstest(i) = u'*Mr*rhs;
%         end
                
        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res; 
    end
        
    if mod(i,5)==0        
        surf(x,y,reshape(Vrp*u,2*N+1,2*N+1))
        view(2)
        shading interp
        title(sprintf('step = %d / %d, rhstest = %g\n',i,Nsteps,rhstest(i)))
        drawnow
    end
end

fprintf('final err = %g\n',sqrt(sum(dx^2*(Vrp*u-Usnap(:,end)).^2)))


