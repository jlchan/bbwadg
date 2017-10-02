function [r s t] = pyr_abctorst(a,b,c)
a = a(:);b = b(:); c = c(:);
r = .5*(1+a).*(1-c) - 1; 
s = .5*(1+b).*(1-c) - 1; 
t = c;
