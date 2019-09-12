d = 3;

N = 1:7;

% 2pt flux ops 
GLL_ops = d*(N+1).^(d+1);
gauss_ops = d*((N+1).^(d+1) + 2*2*(N+1).^d);
gauss_sym_ops = d*((N+1).^(d+1) + 2*(N+1).^d);
stag_ops = d*(N+2).^(d+1);

% print_pgf_coordinates(N,GLL_ops)
% print_pgf_coordinates(N,stag_ops)
% print_pgf_coordinates(N,gauss_ops)
% print_pgf_coordinates(N,gauss_sym_ops)

% matrix ops
mmops = 2*(N+1).^(d+1);
mmops2 = 2*(N+2).^(d+1);
a = 0;
GLL_mops = d*mmops;
gauss_mops = d*(mmops + 2*2*(N+1).^d); % 2x interp, Nfaces_1D, cost
stag_mops = d*mmops2 + 2*2*(N+2).*(N+1).^d; 

% plot(N,GLL_ops,'o-')
% hold on
% plot(N,stag_ops,'s-')
% plot(N,gauss_ops,'x-')

print_pgf_coordinates(N,GLL_mops)
print_pgf_coordinates(N,stag_mops)
print_pgf_coordinates(N,gauss_mops)


plot(N,GLL_mops,'o-')
hold on
plot(N,stag_mops,'s-')
plot(N,gauss_mops,'x-')