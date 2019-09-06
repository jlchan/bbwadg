ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;

N          = 250;
p          = 11;
plotOption = 2;
nPeriods   = 31;
tf         = 100.5;
%closure    = ghostBasis;
closure    = extrapolation;
%closure    = compatibility;

close all
wave( N,p,closure,plotOption,nPeriods,tf );

%%
% plotConvCompatibility
plotConvExtrapolation
% plotConvGhostBasis