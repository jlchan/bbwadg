function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = ...
    euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma,rhoxlog,rhoylog,betaxlog,betaylog)

if nargin == 9
    rholog = logmean(rhox,rhoy);
    betalog = logmean(betax,betay);
else
    rholog = logmean(rhox,rhoy,rhoxlog,rhoylog);
    betalog = logmean(betax,betay,betaxlog,betaylog);
end

% arithmetic avgs
rhoavg = .5*(rhox+rhoy);
uavg = .5*(ux+uy);
vavg = .5*(vx+vy);

vnavg = ux.*uy + vx.*vy; % 2*(uavg.^2 + vavg.^2) - .5*((ux.^2+uy.^2) + (vx.^2+vy.^2));
pa = rhoavg./(betax+betay);
f4aux = rholog./(2*(gamma-1)*betalog) + pa + .5*rholog.*vnavg;

FxS1 = rholog.*uavg;       
FxS2 = FxS1.*uavg + pa;    
FxS3 = FxS1.*vavg;               
FxS4 = f4aux.*uavg;

FyS1 = rholog.*vavg; 
FyS2 = FxS3; %FyS1.*uavg;
FyS3 = FyS1.*vavg + pa;
FyS4 = f4aux.*vavg;

end