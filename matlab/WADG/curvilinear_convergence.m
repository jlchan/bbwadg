% curvilinear p-convergence

L2err1 = [2.429482e-02 4.655834e-03 1.950810e-04 2.093990e-05 8.206020e-07];   % J-dependent mass matrix
L2err2 = [2.429482e-02 4.655300e-03 1.950115e-04 2.093971e-05 8.203261e-07]; % 1/J-weighted
L2err3 = [2.429482e-02 4.655829e-03 1.950574e-04 2.093984e-05 8.205138e-07]; % P_N(1/J) weighted

N = length(L2err1);
semilogy(1:N,L2err1,'o-.','linewidth',2,'markersize',12);
hold on
semilogy(1:N,L2err2,'x--','linewidth',2,'markersize',12)
semilogy(1:N,L2err3,'^:','linewidth',2,'markersize',12)
legend('J-dependent mass matrix','1/J-projection','P_N(1/J)-projection (aliasing)')
set(gca,'fontsize',14)

%% curvilinear h-convergence
rates = [];
N = 2;
err = [ 2.604498e-02  6.615078e-03 7.224455e-04 9.521704e-05]; % J-dependent mass matrix
err2 = [2.604310e-02 6.615064e-03 7.224452e-04 9.521703e-05]; % 1/J weighted
err3 = [2.604313e-02 6.615064e-03 7.224452e-04 9.521703e-05]; % P_N(1/J) weighted
h = .5.^(1:length(err));
loglog(h,err,'o-','linewidth',2,'markersize',12); hold on
loglog(h,err2,'s:','linewidth',2,'markersize',12);
fit = [log(h(:)) ones(size(h(:)))]\log(err(:));
rates = [rates fit(1)];
% print_pgf_coordinates(h,err)
% print_pgf_coordinates(h,err2)
% print_pgf_coordinates(h,err3)
format shorte
fprintf('%4.4e & ',err);fprintf('\n')
fprintf('%4.4e & ',err2);fprintf('\n')
format


N = 3;
err = [1.012391e-02 3.876823e-04 2.888795e-05 1.788749e-06]; % J-dependent mass matrix
err2 = [1.012384e-02 3.876819e-04 2.888795e-05 1.788749e-06]; % 1/J weighted
err3 = [1.012384e-02 3.876819e-04 2.888795e-05 1.788749e-06]; % P_N(1/J) weighted
h = .5.^(1:length(err));
loglog(h,err,'o-','linewidth',2,'markersize',12); hold on
loglog(h,err2,'s:','linewidth',2,'markersize',12);
fit = [log(h(:)) ones(size(h(:)))]\log(err(:));
rates = [rates fit(1)];
% print_pgf_coordinates(h,err)
% print_pgf_coordinates(h,err2)
% print_pgf_coordinates(h,err3)
format shorte
fprintf('%4.4e & ',err);fprintf('\n')
fprintf('%4.4e & ',err2);fprintf('\n')
format

N = 4;
err = [4.831438e-04 4.841638e-05 1.389519e-06 4.715870e-08]; % J-dependent mass matrix
err2 = [4.830959e-04 4.841621e-05 1.389518e-06 4.715869e-08]; % 1/J weighted
err3 = [4.830959e-04 4.841621e-05 1.389518e-06 4.715869e-08]; % P_N(1/J) weighted

h = .5.^(1:length(err));
loglog(h,err,'o-','linewidth',2,'markersize',12); hold on
loglog(h,err2,'s:','linewidth',2,'markersize',12);

fit = [log(h(:)) ones(size(h(:)))]\log(err(:));
rates = [rates fit(1)];
% print_pgf_coordinates(h,err)
% print_pgf_coordinates(h,err2)
% print_pgf_coordinates(h,err3)
format shorte
fprintf('%4.4e & ',err);fprintf('\n')
fprintf('%4.4e & ',err2);fprintf('\n')
format



title(sprintf('Order %d: rate = %f\n',N,fit(1)))
