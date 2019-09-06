N = 1:7;

% num flux evals
d = 3;
tp_GLL = d*(N+1).^(d-1).*(N+1).^2;
tp_GQ = d*(N+1).^(d-1).*((N+1).^2 + 4*(N+1)); % number of lines * 1D nnz
tp_GQ_sym = d*((N+1).^(d-1)).*((N+1).^2 + 2*(N+1));
tp_stag = d*(N+2).^(d-1).*(N+2).^2;

bar(N,[tp_GLL; tp_GQ; tp_GQ_sym; tp_stag]','grouped')
legend('Lobatto','Gauss','Sym Gauss','Staggered','Location','Best')
grid on

return

plot(tp_GLL,'o--')
hold on
plot(tp_GQ,'s--')
plot(tp_GQ_sym,'^--')
plot(tp_stag,'x--')

return

op_GLL = 3*(N+1).^4;
op_GQ = op_GLL + 12*(N+1).^3;
op_stag = 3*(N+2).^4 + 6*(N+1).^3.*(N+2);

plot(op_GLL,'o--')
hold on
plot(op_GQ,'s--')
plot(op_stag,'x--')

% fprintf('\n');
% print_pgf_coordinates(tp_GLL)
% print_pgf_coordinates(tp_GQ)
% print_pgf_coordinates(tp_stag)
fprintf('\n');
print_pgf_coordinates(op_GLL)
print_pgf_coordinates(op_stag)
print_pgf_coordinates(op_GQ)
