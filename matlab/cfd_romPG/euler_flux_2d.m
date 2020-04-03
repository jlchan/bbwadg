function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux_2d(rhoL,rhoR,uL,uR,vL,vR,betaL,betaR)

gamma = 1.4;

% optimized evaluations
rholog = logmean(rhoL,rhoR);
rhoavg = .5*(rhoL+rhoR);
uavg = .5*(uL+uR);
vavg = .5*(vL+vR);
vnavg = uL.*uR + vL.*vR; 
pa = rhoavg./(betaL+betaR);

f4aux = rholog./(2*(gamma-1)*logmean(betaL,betaR)) + pa + .5*rholog.*vnavg;

FxS1 = rholog.*uavg;      
FxS2 = FxS1.*uavg + pa;   
FxS3 = FxS1.*vavg;              
FxS4 = f4aux.*uavg;

FyS1 = rholog.*vavg;
FyS2 = FxS3;
FyS3 = FyS1.*vavg + pa;
FyS4 = f4aux.*vavg;