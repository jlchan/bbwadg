% convergence for spherical solution at t = .5


N = 2;
Np = (N+1)*(N+2)*(N+3)/6;
dofs = [48 385 10087 83762]*Np;
err_skew =   [3.856958e-02 4.487022e-03 7.285771e-05 6.622444e-06];
% err_strong = [3.857076e-02 4.486986e-03 7.285768e-05 6.622444e-06];
err_skew = [1.387606e-02 5.498242e-03 9.033756e-05 7.856624e-06];
% err_strong_underint = [];
loglog(dofs,err_skew,'o-')
hold on
% loglog(dofs,err_strong,'x-')
loglog(dofs,dofs.^(-(N+1)/3),'--')

% print_pgf_coordinates(dofs,err_strong)
h = dofs.^(-1/3); h = h/h(1);
print_pgf_coordinates(h,err_skew)

%%
N = 3;
Np = (N+1)*(N+2)*(N+3)/6;
dofs = [48 385 10087 83762]*Np;
% err_skew =   [6.524365e-03 2.579821e-04 2.194664e-06 1.084056e-07];
% err_strong = [6.521657e-03 2.579401e-04 2.194650e-06 1.084055e-07];
err_skew = [4.076355e-03 3.251153e-04 3.819027e-06 2.078294e-07];

loglog(dofs,err_skew,'o-')
hold on
% loglog(dofs,err_strong,'x-')
loglog(dofs,dofs.^(-(N+1)/3),'--')

print_pgf_coordinates(dofs,err_skew)
% print_pgf_coordinates(dofs,err_strong)

h = dofs.^(-1/3); h = h/h(1);
print_pgf_coordinates(h,err_skew)

%% 
N = 4;
Np = (N+1)*(N+2)*(N+3)/6;
dofs = [48 385 10087]*Np;
% err_skew =  [6.617683e-04 3.599545e-05 8.415544e-08];
err_skew = [8.362368e-04 4.830582e-05 1.007447e-07];
loglog(dofs,err_skew,'o-')
hold on
loglog(dofs,1e3*dofs.^(-(N+1)/3),'--')

print_pgf_coordinates(dofs,err_skew)

h = dofs.^(-1/3); h = h/h(1);
print_pgf_coordinates(h,err_skew)


%% 

N = 5;
Np = (N+1)*(N+2)*(N+3)/6;
dofs = [48 385 10087]*Np;
err_skew            = [1.755100e-04 2.024081e-06 9.371271e-09];
err_strong          = [1.754078e-04 1.973840e-06 9.371272e-09];
err_strong_underint = [1.754320e-04 1.973956e-06 9.371270e-09];
err_skew = [2.327850e-04 3.133785e-06 7.822332e-09];
loglog(dofs,err_skew,'o-')
hold on
loglog(dofs,1e3*dofs.^(-(N+1/2)/3),'--')
polyfit(log(dofs),log(err_skew),1)
%print_pgf_coordinates(dofs,err_skew)
h = dofs.^(-1/3); h = h/h(1);
print_pgf_coordinates(h,err_skew)

