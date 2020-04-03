function [wq id] = get_EQP_nodes(Vtest,wq0,delta,C,d)

if nargin < 4
    C = [];
    d = [];
end

f = ones(size(wq0));
A = [Vtest'
    -Vtest'];

% normalize to get "relative" error
b = [Vtest'*wq0 + delta*max(abs(Vtest'*wq0));
    -Vtest'*wq0 + delta*max(abs(Vtest'*wq0))];

% convergence doesn't improve monotonically as constraint tolerance increases?
%options = optimoptions('linprog','Algorithm','dual-simplex','TolCon',5e-8);
options = optimoptions('linprog','Algorithm','dual-simplex','TolCon',5e-7);
options.Preprocess = 'none'; % from https://www.mathworks.com/matlabcentral/answers/440944-linprog-stopped-because-it-exceeded-its-allocated-memory-with-the-dual-simplex-algorithm

% exact enforcement as Joey recommended
% C = [C;Vtest']; d = [d;Vtest'*wq0]; [wq,~,flag] = linprog(f,[],[],C,d,zeros(size(wq0)),inf(size(wq0)),[],options); % accuracy constraints as hard constraints, lb/ub

[wq,~,flag] = linprog(f,A,b,C,d,zeros(size(wq0)),inf(size(wq0)),[],options); % f, ineq constraints, hard constraints, lb/ub

id = find(abs(wq)>1e-12);
wq = wq(id);

fprintf('%d EQP nodes, residual for linprog = %g, err flag = %d\n',length(wq),norm(Vtest(id,:)'*wq - Vtest'*wq0)/norm(Vtest'*wq0),flag)