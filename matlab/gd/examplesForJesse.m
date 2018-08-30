ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;

N          = 21;
p          = 3;
plotOption = 2;
nPeriods   = 2;
tf         = 1.25;
%closure    = ghostBasis;
closure    = extrapolation;
%closure    = compatibility;

close all
wave( N,p,closure,plotOption,nPeriods,tf );
%%
% plotConvCompatibility
plotConvExtrapolation
% plotConvGhostBasis