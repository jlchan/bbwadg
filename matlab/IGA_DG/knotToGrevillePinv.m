clear
NB = 7;
Ksub = NB;
N = NB+Ksub-1;

VX = linspace(-1,1,Ksub+1);
re = linspace(-1,1,N+1)';
t0 = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]';
t = t0;

load optKnots.mat
if NB > 1
    VXopt = optKnots{NBlist==NB & KsubList==Ksub};
else
    VXopt = linspace(-1,1,Ksub+1);
end

for ii = 1:25
    Ve = bspline_basismatrix(NB+1, t, t0);
    
%     rp = linspace(-1,1,1000)';
%     Vp = bspline_basismatrix(NB+1, t, rp);
    
    t = Ve*re;
    
    if 0
        clf
        plot(rp,Vp)
        hold on
        plot(t0,Ve,'o','markersize',8,'linewidth',2)
        plot(t0,(1+t)/2,'x--','markersize',8)
        plot(t,-.025 + t*0,'o')
        
        Vpopt = bspline_basismatrix(NB+1, topt, rp);
        plot(rp,Vpopt,'--')
        plot(t0,(1+topt)/2,'*--','markersize',8)
        plot(topt,-.025 + topt*0,'x')
        
        pause
    end
end

% knot to greville operator
G = zeros(NB+Ksub,2*NB+Ksub+1);
for i = 1:NB+Ksub
    G(i,i+1:i+NB) = 1 / NB;
end

topt = [VXopt(1)*ones(1,NB) VXopt VXopt(end)*ones(1,NB)]';
r0 = G*t0;
r = G*t;
ropt = G*topt;

g0 = sum(G(:,end-NB:end),2)-sum(G(:,1:NB+1),2); % contribution from fixed knot locations
Gred = G(:,NB+2:end-NB-1);
ids = NB+2:NB+Ksub;
norm(ropt - (Gred*topt(ids) + g0))
norm(topt(ids) - pinv(Gred)*(ropt-g0))

% t = t0; 
% for ii = 1:25
%     Ve = bspline_basismatrix(NB+1, t, r0); % eval at original greville pts
%     rt = Ve*re; % smooth greville points - assume abscissae are linearly spaced
%     
%     a = .33;
%     rt = a*rt + (1-a)*r0;
%     
%     % recover knots from greville
%     t(ids) = pinv(Gred)*(rt-g0);
%     
%     t = sort(t,'ascend');    
%     
%     clf
%     plot(topt,topt*0,'o');
%     hold on
%     plot(t,t*0,'x');
%     pause
% end
% rt = G*t;

hold on;
plot(t,t*0,'bo','markersize',8)
plot(topt,topt*0,'rx','markersize',8);
% plot(rt,rt*0,'kd','markersize',8)
