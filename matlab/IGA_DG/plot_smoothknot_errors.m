% recompute spline knots with new eigenvalue problem

% N = 2
KK = [4 8 16 32 64];
terr = [0.0121 0.0075 0.0046 0.0023 0.0012];
rerr = [0.0060 0.0038 0.0024 0.0011 5.8806e-04];
plot(KK,rerr,'o--')
print_pgf_coordinates(KK,rerr)
hold on

% N = 3
terr = [0.0090 0.0041 0.0028 0.0020 0.0012];
rerr = [0.0030 0.0020 0.0018 0.0012 6.7857e-04];
print_pgf_coordinates(KK,rerr)
plot(KK,rerr,'o--')
 
% N = 4
terr = [0.0016  0.0101 0.0074 0.0048 0.0028];
rerr = [4.01e-4 0.0025 0.0016 8.9352e-04 6.3901e-04];
print_pgf_coordinates(KK,rerr)
plot(KK,rerr,'o--')


% % N = 5
% terr = [0.0276 ];
% rerr = [0.0055 ];
% plot(KK,terr,'o--')

legend('N = 2','N = 3','N = 4','N = 5')