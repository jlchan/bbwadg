function [r s t] = wedge_abctorst(a,b,c)

% transform to wedge 
r = .5*(1+a).*(1-c) - 1;
s = b; 
t = c;