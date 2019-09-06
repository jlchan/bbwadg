% load Usnap_euler
load Usnap_euler_wall

gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

K = size(Usnap(:,1),1)/3;
U0 = reshape(Usnap(:,1),K,3);
%  rho = Vr'*rho;
%  m = Vr'*m;
%  E = Vr'*E;

Nsample = 10; %ceil(size(Usnap,2)/500); % downsample
Us1 = Usnap((1:K),1:Nsample:end);
Us2 = Usnap((1:K)+K,1:Nsample:end);
Us3 = Usnap((1:K)+2*K,1:Nsample:end);

% add entropy variables to snapshots
Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
[VrUV,SrUV,~] = svd(Us,0);

Us = [Us1 Us2 Us3];
[VrU,SrU,~] = svd(Us,0);

%% component singular value decay

figure
semilogy(diag(SrU),'o','linewidth',2,'markersize',10)
hold on
semilogy(diag(SrUV),'.-','linewidth',2,'markersize',10)

h = legend('Without entropy variables','With entropy variables','Location','Best');
set(h,'fontsize',16)
set(gca,'fontsize',15)
grid on
xlabel('Mode index','fontsize',15)
xlim([0,900])

%% vector vars instead of components

Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
[VrUc,SrUc,~] = svd(Us,0);

UU = [Us1; Us2; Us3];
VV = [V1(Us1,Us2,Us3); V2(Us1,Us2,Us3); V3(Us1,Us2,Us3)];
[VrUvec,SrUvec,~] = svd(UU,0);
[VrVvec,SrVvec,~] = svd(VV,0);
[VrUVvec,SrUVvec,~] = svd([UU VV],0);

sUc = diag(SrUc);
sU = diag(SrUvec);
sV = diag(SrVvec);
sUV = diag(SrUVvec);

sUc_energy = sUc./sum(sUc);
sU_energy = sU./sum(sU);
sV_energy = sV./sum(sV);
sUV_energy = sUV./sum(sUV);

figure
semilogy(sU_energy,'o-','linewidth',2,'markersize',10)
hold on
semilogy(sV_energy,'.--','linewidth',2,'markersize',10)
semilogy(sUV_energy,'x--','linewidth',2,'markersize',10)
semilogy(sUc_energy,':','linewidth',2,'markersize',10)
legend('U only','V only','U and V','U+V components')
set(gca,'fontsize',15)
grid on
xlabel('Mode index','fontsize',15)


%%

Nmodes = 25;
PU = VrU(:,1:Nmodes)*VrU(:,1:Nmodes)';
PUV = VrUV(:,1:Nmodes)*VrUV(:,1:Nmodes)';
for i = 1:size(Us1,2)
    f1 = Us1(:,i);
    f2 = Us2(:,i);
    f3 = Us3(:,i);
    errU(i) = sqrt(norm(f1-PU*f1)^2+norm(f2-PU*f2)^2+norm(f3-PU*f3)^2);
    errUV(i) = sqrt(norm(f1-PUV*f1)^2+norm(f2-PUV*f2)^2+norm(f3-PUV*f3)^2);
end
%%
sk = 5;
semilogy(1:sk:length(errUV),errUV(1:sk:end),'o-','linewidth',2,'markersize',10)
hold on
semilogy(1:sk:length(errUV),errUV(1:sk:end),'x--','linewidth',2,'markersize',10)
semilogy(1:sk:length(errUV),abs(errU(1:sk:end)-errUV(1:sk:end)),'*--','linewidth',2,'markersize',10)
set(gca,'fontsize',15)
grid on
xlabel('Snapshot index','fontsize',15)
h = legend('Without entropy variables','With entropy variables','Difference in errors','Location','Best');
set(h,'fontsize',16)

%%
figure
x = linspace(-1,1,K)';
sk = 50;
for id = 5;%1:size(VrU,2)
    clf
    [~,maxid] = max(abs(VrU(:,id)));
    scale = VrUV(maxid,id)/VrU(maxid,id);
    plot(x(1:sk:end),scale*VrU((1:sk:end),id)/max(abs(VrU(:,id))),'o-','linewidth',2,'markersize',10)
    hold on
    plot(x(1:sk:end),VrUV((1:sk:end),id)/max(abs(VrUV(:,id))),'x--','linewidth',2,'markersize',10)
%     pause
end
h = legend('Without entropy variables','With entropy variables','Location','Best');
set(h,'fontsize',16)
set(gca,'fontsize',15)
grid on
% xlabel('Snapshot index','fontsize',15)
