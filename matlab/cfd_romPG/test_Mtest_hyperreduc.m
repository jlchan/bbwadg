K = 2500;

Nmodes = 75;

xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;
S = sparse(S);

load Usnap_euler_wall.mat

U0 = reshape(Usnap(:,1),K,3);

Nsample = 1; %ceil(size(Usnap,2)/500); % downsample
Us1 = Usnap((1:K),1:Nsample:end);
Us2 = Usnap((1:K)+K,1:Nsample:end);
Us3 = Usnap((1:K)+2*K,1:Nsample:end);

gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

% add entropy variables to snapshots
Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
tic;
[Vr,Sr,~] = rsvd(Us,5*Nmodes);
fprintf('Time for generating Vr = %f\n',toc)

Vr = Vr(:,1:Nmodes);
Vrfull = Vr;

sig = diag(Sr);
tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2));

% add range + constants to snapshots
tic;[Vtest Stest, ~] = svd([ones(size(x)) Vr S*Vr],0);
fprintf('Time for generating Vtest = %f\n',toc)
sigt = diag(Stest);
Vtest = Vtest(:,sigt > 1e-12);

%%

% empirical cubature
Va = Vr;
Vb = Vr; 
tic;
Vmass = zeros(size(Va,1),size(Va,2)*(size(Vb,2)+1)/2);
sk = 1;
for i = 1:size(Va,2)
    for j = i:size(Vb,2)
        Vmass(:,sk) = Va(:,i).*Vb(:,j);
        sk = sk + 1;
    end
end
fprintf('Time for generating Vmass = %f\n',toc)

% reduce target space
tic; [Vmass,Smass,~] = rsvd(Vmass,25*Nmodes); toc
smass = diag(Smass);
smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
Vmass = Vmass(:,smass_energy > tol);

%% compute vol points

w0 = ones(K,1)*dx;
[wr id] = get_empirical_cubature(Vmass,w0,tol,25*Nmodes);

% make new mass, projection, interp matrices
Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));

%% compute stabilizing points for gappy POD

sk = 1;
while cond(Vtest(id,:)) > 1e10
    ids = setdiff(1:K,id);
    for i = 1:length(ids)
        idi = [id; i];
        condV(i) = cond(Vtest(idi,:)'*Vtest(idi,:));
        if mod(i,100)==0
            fprintf('on i = %d out of %d\n',i,length(ids))
        end
    end
    [~,i] = min(condV);
    id = [id;ids(i)];
    sk = sk + 1;
end

return

%% compute stabilizing points for quadrature

Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
[Vx Dx] = eig(Mtest);
d = diag(Dx);
[d,p] = sort(d,'ascend');
Z = Vtest*Vx(:,p(d < 1e-12)); 
    
Zmass = zeros(size(Z,1),size(Z,2)*(size(Z,2)+1)/2);
sk = 1;
for i = 1:size(Z,2)
    for j = 1:size(Z,2)
        Zmass(:,sk) = Z(:,i).*Z(:,j); 
        sk = sk + 1;
    end
end
[Zmass,SZ,~] = svd(Zmass,0); 
Zmass = Zmass(:,diag(SZ) > 1e-12);
fprintf('\n getting additional points for stability\n')
[wz idz] = get_empirical_cubature(Zmass,w0,tol,size(Zmass,2));

idt = unique([id; idz]);

dz = .001;
for zscale = dz:dz:.1
    J = [Vmass zscale*Zmass];
    b = sum(J',2);
    w_unscaled = lsqlin(J(idt,:)',b,[],[],[] ,[] ,zeros(size(idt)),inf(size(idt)),J(idt,:)'\b); % upper bound = total vol
    w = w_unscaled*sum(w0)/K;
    
    Mtest = Vtest(idt,:)'*diag(w)*Vtest(idt,:);
    [Vt Dt] = eig(Mtest);
    lam = diag(Dt);
    [~,idx] = min(lam);
    
    e1 = norm(Vr(id,:)'*diag(wr)*Vr(id,:) - dx*eye(Nmodes),'fro');
    et = norm(Vr(idt,:)'*diag(w)*Vr(idt,:) - dx*eye(Nmodes),'fro');
    fprintf('\n original HR err = %g, new HR err = %g. min eig/dx = %g\n',e1,et,min(lam)/dx)
    
    semilogy(zscale,et,'bo')
    hold on
    semilogy(zscale,min(lam)/dx,'rx')    
end

return
%%

% plot(x(id),wr,'o')
% hold on
% plot(x(idt),w,'x')
% return
% clf
plot(x,Vtest*Vt(:,idx),'o-')
hold on
plot(x(idz),Vtest(idz,:)*Vt(:,idx),'x','linewidth',2)
