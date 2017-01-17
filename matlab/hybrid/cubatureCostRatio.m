N = 3;
NpH = (N+1).^3;
NpW = (N+1).^2.*(N+2)/2;
NpP = (N+1).*(N+2).*(2*N+3)/6;
NpT = (N+1).*(N+2).*(N+3)/6;

for NN = N
    r = tet_cubature(2*NN);
    NcT(NN) = length(r);
end

plot(NcT./NpT,'.-')
hold on
Np3 = (N+1).^3;
plot(Np3./NpW,'r.-')
plot(Np3./NpP,'k.-')
