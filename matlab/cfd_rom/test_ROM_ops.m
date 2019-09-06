K = 1000;
FinalTime = .3;
xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

e = ones(K-1,1);
Q = diag(e,1)-diag(e,-1);
Q(1,:) = 0; Q(:,1) = 0;
Q(end,:) = 0; Q(:,end) = 0;
Q(1,end) = -1; Q(end,1) = 1;
Q(1,2) = 1; Q(2,1) = -1;
Q(K-1,K) = 1; Q(K,K-1) = -1;

Q = sparse(Q);
KS = sparse(2*eye(K) - abs(Q))/dx;

%% entropy vars 

gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

%% construct ROM

load Usnap_euler
% load Fsnap_euler

U0 = reshape(Usnap(:,1),K,3);
%  rho = Vr'*rho;
%  m = Vr'*m;
%  E = Vr'*E;

Nsample = 2; %ceil(size(Usnap,2)/500); % downsample
Us1 = Usnap((1:K),1:Nsample:end);
Us2 = Usnap((1:K)+K,1:Nsample:end);
Us3 = Usnap((1:K)+2*K,1:Nsample:end);

% add entropy variables to snapshots
Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
%     Us = [Us1 Us2 Us3];
[Vr,Sr,~] = svd(Us,0);

Nmodes = 25;
Vr = Vr(:,1:Nmodes);
Vrp = Vr;

sig = diag(Sr);
tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))

% add range + constants to snapshots
[Vtest stest, ~] = svd([Vr Q*Vr]);
sigt = diag(stest);
Vtest = Vtest(:,sigt > 1e-10);
Vtest = orth([ones(size(x)) Vtest]); % ensure 1
[Vrange, Srange, ~] = svd(Q*Vr,0);
Vrange = Vrange(:,diag(Srange)>1e-10);

Qfull = Q;
Kfull = KS;

rho = Vr'*U0(:,1);
m = Vr'*U0(:,2);
E = Vr'*U0(:,3);

%     % set tol based on init condition
%     tol = sqrt(sum(sum(dx*(U0-Vrp*[rho m E]).^2)))/sqrt(sum(sum(dx*U0.^2)))

%% empirical cubature

% target space
Vtest1 = Vr;
Vtest2 = Vr;
% Vtest2 = Vtest;
sk = 1;
Vmass = zeros(size(Vtest1,1),size(Vtest1,2)*size(Vtest2,2));
for i = 1:size(Vtest1,2)
    for j = 1:size(Vtest2,2)
        Vmass(:,sk) = Vtest1(:,i).*Vtest2(:,j);
        sk = sk + 1;
    end
end

% reduce target space
[Vmass,Qmass,~] = svd(Vmass,0);
smass = diag(Qmass);
Qmass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
% Vmass = Vmass(:,find(Qmass_energy > tol));
% Vmass = Vmass(:,1:size(Vtest,2)*2);

% semilogy(smass,'bo')
% hold on
% id = find(Qmass_energy > tol);
% semilogy(id,smass(id),'rs','linewidth',2,'markersize',12)
% id = size(Vtest,2)*2;
% semilogy(id,smass(id),'k*','linewidth',2,'markersize',12)
% return

Nred = find(Qmass_energy < tol,1,'first');
% Nred = find(Qmass_energy < tol^2,1,'first');
Vmass = Vmass(:,1:max(Nred,size(Vtest,2)*2));
% Vmass = Vmass(:,1:Nred);
% Vmass = Vmass(:,1:size(Vtest,2)*2);   

[wr id] = get_empirical_cubature(Vmass,ones(size(x(:)))*dx,tol,inf); % cheaper
% wr = ones(size(x(:))); id = 1:length(x(:));
%  [wr id] = get_EQP_nodes(Vmass,ones(size(x(:)))*dx,tol);

% norm(dx*sum(Vmass',2) - Vmass(id,:)'*wr)
% norm(dx*sum(Vmass',2) - Vmass(id,:)'*wr,'inf')

% make new mass, projection, interp matrices
Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
Pr = Mr\(Vr(id,:)'*diag(wr));
% Vr = Vr(id,:);

Ptest = (Vtest(id,:)'*diag(wr)*Vtest(id,:))\(Vtest(id,:)'*diag(wr));
% Vr1 = [Vr ones(size(Vr,1),1)];
% Pr1 = (Vr1(id,:)'*diag(wr)*Vr1(id,:))\(Vr1(id,:)'*diag(wr));
% R = Vtest\Vr1;
% Ptest = R*Pr1;
Q = Ptest'*(Vtest'*Qfull*Vtest)*Ptest;
KS = Ptest'*(Vtest'*Kfull*Vtest)*Ptest;

PtestFull = zeros(size(Vtest'));
PtestFull(:,id) = Ptest;
Qs = PtestFull'*(Vtest'*Qfull*Vtest)*PtestFull;

%%
% approx err
% u = Vtest*(1:size(Vtest,2))';
% v = Vr*(1:size(Vr,2))';

u = Vr(:,12); u = u/max(abs(u));
v = Vr(:,13); v = v/max(abs(v));
norm(Vr'*Qfull*u - Vr(id,:)'*Q*u(id))

% norm(Qfull*u - Qs*u)
norm(Qfull*u-Qs*u)

norm(Vr'*(v.*(Qfull*u)) - Vr(id,:)'*(v(id).*(Q*u(id)))) % around the residual of emp. cubature




