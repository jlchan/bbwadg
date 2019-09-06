clear
K = 250;
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
Kr = (2*eye(size(S))-abs(S))/dx;
Krfull = Kr;
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
% fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);
fD = @(uL,uR) .5*(abs(uL)+abs(uR)).*(uR-uL);
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

u0 = .0-sin(pi*x);

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
    u = Pr*u0;
    
else
    %     load Usnap_burgers_moving
    load Usnap_burgers_GLF
    %     load Usnap_burgers_LLF
    %     load Usnap_burgers
    [Ur,Sr,Vr] = svd(Usnap,0);
    
    Nmodes = 3;
    Vr = [x(:).^0 Ur(:,1:Nmodes)];
    
    Vrp = Vr;
    
    invM = inv(dx*Vr'*Vr);
    Pr = invM*Vr'*dx;
    
    % project solution
    u = Pr*u0;
    
end



opt = 1;
if opt==0
    % no hyperreduction
    id = 1:size(Vrp,1);
    idD = id;
    
    Sfull = S;
    [Utest,Stest,~] = svd([ones(K,1) Vr S*Vr],0); % enrich with "least squares" basis
    Vtest = Utest(:,diag(Stest) > 1e-11); % include nonzero energy modes
    S = Vtest*Vtest'*Sfull*Vtest*Vtest'; % project matrix
    
    % viscosity matrix
    Kr = Vr'*Kr*Vr;
    
elseif opt==1
    
    Sfull = S;
    [Utest,Stest,~] = svd([ones(K,1) Vr S*Vr],0); % enrich with "least squares" basis
    %     [Utest,Stest,~] = svd([ones(K,1) S*Vr],0); % enrich with "least squares" basis
    Vtest = Utest(:,diag(Stest) > 1e-11); % include nonzero energy modes
    
    [Vrange,Srange,~] = svd(S*Vr,0);
    Vrange = Vrange(:,diag(Srange)>1e-11);
    sig = diag(Stest);
    tol = sum(sig(Nmodes+1:end).^2)/sum(sig.^2)
    %     tol = sqrt(sum(dx*(Vr*u - u0).^2))/sqrt(sum(dx*u0.^2))
    
    V1 = Vr;    
    V2 = Vtest;
%     V2 = Vr;
    sk = 1;
    for i = 1:size(V1,2)
        for j = 1:size(V2,2)
            Vmass(:,sk) = V1(:,i).*V2(:,j);
            sk = sk + 1;
        end
    end
    
    % reduce target space
    [Vmass,Smass,~] = svd(Vmass,0);
    smass = diag(Smass);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    %Vmass = Vmass(:,smass_energy > tol);    
    Vmass = Vmass(:,smass_energy > 1e-13);    
    disp('Done with Vmass svd')
    
    [wr id] = get_empirical_cubature(Vmass,ones(K,1)*dx,tol,25*Nmodes);
%     [wr id] = get_empirical_cubature(Vmass,ones(K,1)*dx,1e-13,25*Nmodes);
    %     wr = ones(K,1)*dx;
    %     id = 1:K;
        
%     [Vr1, S1,~] = svd([Vr ones(size(Vr,1),1)],0); 
%     Vr1 = Vr1(:,diag(S1)>1e-13);
%     Pr1 = (Vr1(id,:)'*diag(wr)*Vr1(id,:))\(Vr1(id,:)'*diag(wr));
%     R = Vtest\Vr1;
%     Ptest = R*Pr1;
    Ptest = (Vtest(id,:)'*diag(wr)*Vtest(id,:))\(Vtest(id,:)'*diag(wr));
    Stest = Ptest'*(Vtest'*Sfull*Vtest)*Ptest;
    
    keyboard
    
    % viscosity matrix
    Kr = Vr'*Kr*Vr;
    
    %     idD = id;
    %     idD = uniquetol([id-1,id,id+1]);
    %     idD(idD < 1 | idD > K) = [];
    
    %     SD = abs(Sfull(idD,idD));
    
    %     VrD = Vr(idD,:);
    Vrfull = Vr;
    Vr = Vr(id,:);
    
    % hyperreduction
    S = Stest;
    %     SD = abs(Sfull);
    
    idD = id;
    idD = unique([id-1;id;id+1]); idD = max(min(idD,K),1);
    idD = 1:size(Sfull,1);
    SD = abs(Sfull(idD,idD));
    
    %     wrD = lsqlin(Vmass(idD,:)',sum(Vmass'*dx,2),[],[],[],[],zeros(size(idD)),[],Vmass(idD,:)'\(sum(Vmass'*dx,2))); % upper bound = total vol
    %     Prp = (Vtest(idD,:)'*diag(wrD)*Vtest(idD,:))\(Vtest(idD,:)'*diag(wrD));
    %     SD = Prp'*Vtest'*abs(Sfull)*Vtest*Prp;
    
end

dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));

for i = 1:Nsteps
    for INTRK = 1:5
        [ux uy] = meshgrid(Vr*u);
        rhs = Vr'*(sum(S.*fS(ux,uy),2));
        
        if (INTRK==5)
            rhstest(i) = u'*rhs;
        end
        
        tau = 5e-3;
        rhs = (rhs + tau*Kr*u); % invert mass + add art visc
        rhs = -invM*rhs;
        
        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res;
    end
    
    if mod(i,5) == 0
        %         plot(x,sum(SD.*fD(ux,uy),2),'o--')
        %         hold on
        %         plot(x,Vtest*Vtest'*sum(SD.*fD(ux,uy),2),'x--')
        %         plot(x,dx*Krfull*Vrp*u,'ks--')
        plot(x(id),Vr*u,'bo')
        hold on
        plot(x,Vrp*u,'b.')
        plot(x,Usnap(:,i),'r-')
        hold off
        %         axis([-1,1,-.1,.1])
        title(sprintf('max rhstest = %g\n',max(abs(rhstest(i)))))
        
        drawnow
    end
end

err = sqrt(sum(sum(dx*(Vrp*u - Usnap(:,i)).^2)))
title(sprintf('Error = %g\n',err))
return

%%
plot(x(id),Vr*u,'bo')
hold on
plot(x,Vrp*u,'b.')
plot(x,Usnap(:,end),'r-')
hold off
axis([-1,1,-1,2.5])
title(sprintf('L2 error = %g\n',sqrt(sum(dx*abs(Vrp*u-Usnap(:,end)).^2))))
