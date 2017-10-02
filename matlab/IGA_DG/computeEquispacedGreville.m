function computeEquispacedGreville

Ksub = 2;
NB = 2;
N = NB+Ksub-1;

VX = linspace(-1,1,Ksub+1);
VX0 = VX;
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
r = greville(NB,VX);
% keyboard

% smooth knots
for i = 1:1
    re = linspace(-1,1,N+1)';
    reKsub = linspace(-1,1,Ksub+1)';
    Ve = bspline_basismatrix(NB+1, t, reKsub);
    VX = Ve*re; VX = VX(:)';
    
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        rnew(i,1) = mean(t((i+1):(i+NB))); % greville
    end
end

makeVX = @(VXint) [-1 VXint 1];
VXint = fsolve(@(VXint) sum((greville(NB,makeVX(VXint)) - linspace(-1,1,N+1)').^2),VX(2:end-1));
VX2 = makeVX(VXint);

keyboard
return

% rp = linspace(-1,1,500);
% Vp = bspline_basismatrix(NB+1, t, rp);
% plot(rp,Vp*re,'--');
% hold on
% plot(VX,VX*0,'o');
% hold on
% plot(VX0,VX0*0,'d');return
%
% figure
plot(r,r*0,'o','markersize',8)
hold on
plot(rnew,rnew*0,'d','markersize',8)
% legend('Original Greville','New Greville')

% figure
% plot(VX0,VX0-VX,'o--') % plot knot displacement
hold on

function r = greville(NB,VX)

Ksub = length(VX)-1;
N = NB + Ksub - 1;
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
for i = 1:N+1
    r(i,1) = mean(t((i+1):(i+NB))); % greville
end