N = 1:7;

tp_GLL = 3*(N+1).^4;
tp_GQ = tp_GLL + 12*(N+1).^3;
tp_stag = 3*(N+2).^4;

% plot(tp_GLL,'o--')
% hold on
% plot(tp_GQ,'s--')
% plot(tp_stag,'x--')

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
