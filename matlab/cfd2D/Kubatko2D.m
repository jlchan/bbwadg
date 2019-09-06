function [rq,sq,wq] = Kubatko2D( N, opt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load QuadratureRules.mat;
quadpts = Q_GaussLegendre{N};
if strcmp(opt,'Lobatto')==1
    quadpts = Q_GaussLobatto{N};
end
rs = quadpts.Points;
wq = quadpts.Weights;
rq = rs(:,1);
sq = rs(:,2);

end

