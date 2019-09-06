% SBP spectral method

clear
N = 50; % must be even

Np1D = N;

dx = 2/N;
h = 2*pi/N; % original h on [-pi,pi]
x = linspace(-1,1,N+1); x = x(1:end-1)';
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column,column([1 N:-1:2]))*pi;
W = dx*eye(Np1D);

Q = W*D;
K = D'*W*D;

%% make 2D

[x y] = meshgrid(x);
Np = Np1D^2;


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

Nmodes = 18;
tol = sqrt(sum(sig(Nmodes+1:end).^2)./sum(sig.^2));
tol = max(tol,5e-7); 
 
Vr = Vs(:,1:Nmodes);
Vrp = Vr;

% full ops
Qxfull = kron(Q,W);
Qyfull = kron(W,Q);
Kfull = kron(K,W) + kron(W,K); % d2u/dx2 + d2u/dy2

% test space
[Vtest, Stest, ~] = svd([Vr Qxfull*Vr Qyfull*Vr],0);
% [Vtest, Stest, ~] = svd([Vr (Qxfull+Qyfull)*Vr],0);
sigt = diag(Stest);
sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
Vtest = orth([ones(Np,1) Vtest(:,sigerr > 1e-12)]);

[Vtestx, Stest, ~] = svd([Vr Qxfull*Vr],0);
sigt = diag(Stest);
sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
Vtestx = orth([ones(Np,1) Vtestx(:,sigerr > 1e-12)]);

[Vtesty, Stest, ~] = svd([Vr Qyfull*Vr],0);
sigt = diag(Stest);
sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
Vtesty = orth([ones(Np,1) Vtesty(:,sigerr > 1e-12)]);

u = Vr'*U0(:,1);
% fprintf('tol = %g, initial err = %g\n',tol, sqrt(sum(dx^2*(Vrp*u-U0).^2)))


%% hyperreduc

if 1
    
%     [Vtest,Stest,~] = svd([Vtestx Vtesty],0);
%     stest = diag(Stest);
%     stest_energy = sqrt(1 - (cumsum(stest.^2)./sum(stest.^2)));
%     %Vtest = Vtest(:,stest_energy > tol);
%     Vtest = Vtest(:,stest_energy > 1e-8);
    
    clear Vmass
    sk = 1;
    for i = 1:size(Vr,2)
        for j = 1:size(Vtest,2)
            Vmass(:,sk) = Vr(:,i).*Vtest(:,j);
            sk = sk + 1;
        end
    end
    
    [Vmass,Smass,~] = svd(Vmass,0);
    smass = diag(Smass);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    Vmass = Vmass(:,smass_energy > 1e-8);
%     Vmass(:,1:8*Nmodes);
%     keyboard
    
    maxpts = 10*Nmodes;
    [wr id] = get_empirical_cubature(Vmass,ones(size(x(:)))*dx^2,tol,maxpts);        
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));    
    Pr = Mr\(Vr(id,:)'*diag(wr));
    invMVrT = Mr\(Vr(id,:)');
    Vr = Vr(id,:);
    
    % make new ops
%     Ptestx = (Vtestx(id,:)'*diag(wr)*Vtestx(id,:))\(Vtestx(id,:)'*diag(wr));        
%     Ptesty = (Vtesty(id,:)'*diag(wr)*Vtesty(id,:))\(Vtesty(id,:)'*diag(wr));    
%     Qx = Ptestx'*(Vtestx'*Qxfull*Vtestx)*Ptestx;
%     Qy = Ptesty'*(Vtesty'*Qyfull*Vtesty)*Ptesty;
    Ptest = (Vtest(id,:)'*diag(wr)*Vtest(id,:))\(Vtest(id,:)'*diag(wr));        
    Qx = Ptest'*(Vtest'*Qxfull*Vtest)*Ptest;
    Qy = Ptest'*(Vtest'*Qyfull*Vtest)*Ptest;
    K = Pr'*(Vrp'*Kfull*Vrp)*Pr;
    

    %     K = Vr'*Kfull*Vr;
else
    [Vtest, Stest, ~] = svd([Vr Qxfull*Vr Qyfull*Vr],0);
    sigt = diag(Stest);
    sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
    Vtest = orth([ones(Np,1) Vtest(:,sigerr > 1e-12)]);

    % % full ops
    Qx = Vtest*(Vtest'*Qxfull*Vtest)*Vtest';
    Qy = Vtest*(Vtest'*Qyfull*Vtest)*Vtest';
    K = Vtest*(Vtest'*Kfull*Vtest)*Vtest';
    invMVrT = (1/dx^2)*Vr'; % orth basis
    Pr = Vr';
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
        
        tau = .25*dx;
        rhs = -invMVrT*(rhs + tau*K*uq);
        
%         if INTRK==5
%             rhstest(i) = u'*Mr*rhs;
%         end
                
        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res; 
    end
        
    if mod(i,5)==0        
        surf(x,y,reshape(Vrp*u,Np1D,Np1D))
        view(2)
        shading interp
        title(sprintf('step = %d / %d, rhstest = %g\n',i,Nsteps,rhstest(i)))
        drawnow
    end
end

err = Vrp*u-Usnap(:,end);
figure
surf(x,y,reshape(err,Np1D,Np1D));shading interp
fprintf('final err = %g\n',sqrt(sum(dx^2*(err).^2)))


