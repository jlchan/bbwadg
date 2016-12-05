B1 = [1 0; 0 0; 0 1];
B2 = [0 0; 0 1; 1 0];

A1 = [zeros(2) B1';
    B1 zeros(3)];
A2 = [zeros(2) B2';
    B2 zeros(3)];

n = randn(2,1); n = n/norm(n);
% n = [-1;0];
An = A1*n(1) + A2*n(2);
eig(An)
sqrt(1 + abs(n(1)*n(2)))

%% 3d 

B1 = zeros(6,3); B1(1,1) = 1; B1(6,2) = 1; B1(5,3) = 1;
B2 = zeros(6,3); B2(6,1) = 1; B2(2,2) = 1; B2(4,3) = 1;
B3 = zeros(6,3); B3(5,1) = 1; B3(4,2) = 1; B3(3,3) = 1;

A1 = [zeros(3) B1';
    B1 zeros(6)];
A2 = [zeros(3) B2';
    B2 zeros(6)];
A3 = [zeros(3) B3';
    B3 zeros(6)];

n = randn(3,1); 
% n = [1;1;1];
n = n/norm(n);
syms nx ny nz
n = [nx; ny; nz];

% An = [1            n1*n2,              n1*n3]
%  [    n1*n2,       1              n2*n3,       
%  [    n1*n3,         n2*n3,         1 
An = A1*n(1) + A2*n(2) + A3*n(3);

return
E = An*An; E = E(1:3,1:3); E = E - diag(diag(E));

norm(An)
sqrt(norm(E)+1)
% 1+sum(abs(E),2)
% sqrt(1 + 2*n(1)*n(2))
% sqrt(1 + norm(n,1))
