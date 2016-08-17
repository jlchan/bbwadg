function [r s w a b] = tri_tp_cubature(Nq)

[a,wa] = JacobiGQ(0,0,Nq);
[b,wb] = JacobiGQ(0,0,Nq);
a = ones(length(a),1)*a';
b = b*ones(1,length(b));

r = 0.5*(1+a).*(1-b)-1;
s = b;
w = 0.5*wb*(wa');
w = w.*(1-b); % if using GQ all around

% r = a;
% s = 0.5*(1+b).*(1-a)-1;
% w = 0.5*wb*(wa');
% w = w.*(1-a); % if using GQ all around


r = r(:); 
s = s(:);
w = w(:);