function [H,Q,X]=CSBPp4(n)
%n = number of nodes in the nodal distribution

Q = zeros(n,n);
H = eye(n,n);
X = zeros(n,1); 

%construct the nodal distribution
dx = 2/(n-1);

for i = 1:n
  X(i) = -1+(i-1)*dx; 
end

Hdiag =  [0.294890676177878558830939783321e0 0.152572062389770723104056437390e1 0.257452876984126984126984126984e0 0.179811370149911816578483245150e1 0.412708057760141093474426807760e0 0.127848462301587301587301587302e1 0.923295579805996472663139329806e0 0.100933386085915847820609725372e1];

% construct the H matrix_type
for i = 1:8
H(i,i)=Hdiag(i);
H(n-i+1,n-i+1) = Hdiag(i);
end

% add the mesh spacing
H = dx*H;

% construct the upper portion of Q
Q(1:8,1:12) =    [-0.500000000000000000000000000000e0 0.666586943328533355994182757055e0 -0.244487517322560325374506836830e-1 -0.215737532995325199481303844018e0 0.168150575610263709309742692500e-1 0.728739975272214371555571457780e-1 -0.541161324401972128760463047040e-2 -0.106781004451802107743550139118e-1 0 0 0 0; -0.666586943328533355994182757055e0 0 0.177599254978796912029901169234e0 0.701027851062354857318848322706e0 -0.667819656104436059449690976518e-1 -0.174951380868956741222460493908e0 -0.147104117213342822635593073689e-2 0.311642249389153620392187874115e-1 0 0 0 0; 0.244487517322560325374506836830e-1 -0.177599254978796912029901169234e0 0 0.156293220071511060156262128132e0 0.543917056186048188341586152164e-1 -0.108907345408875248849424004693e0 0.664680844387540608615917538523e-1 -0.150951614734538115101380069557e-1 0 0 0 0; 0.215737532995325199481303844018e0 -0.701027851062354857318848322706e0 -0.156293220071511060156262128132e0 0 0.240849527296313714359679238313e0 0.531915236612749548912014102342e0 -0.107022280849150358362624702355e0 -0.241589449213721869152620314805e-1 0 0 0 0; -0.168150575610263709309742692500e-1 0.667819656104436059449690976518e-1 -0.543917056186048188341586152164e-1 -0.240849527296313714359679238313e0 0 0.312202228759906012618089180688e0 -0.883136973085797776826612695337e-1 0.249572219856036346729865425447e-1 -0.357142857142857142857142857143e-2 0 0 0; -0.728739975272214371555571457780e-1 0.174951380868956741222460493908e0 0.108907345408875248849424004693e0 -0.531915236612749548912014102342e0 -0.312202228759906012618089180688e0 0 0.732608361099657538517277492903e0 -0.133999434001422053713025372219e0 0.380952380952380952380952380952e-1 -0.357142857142857142857142857143e-2 0 0; 0.541161324401972128760463047040e-2 0.147104117213342822635593073689e-2 -0.664680844387540608615917538523e-1 0.107022280849150358362624702355e0 0.883136973085797776826612695337e-1 -0.732608361099657538517277492903e0 0 0.762334003440718790010098904135e0 -0.200000000000000000000000000000e0 0.380952380952380952380952380952e-1 -0.357142857142857142857142857143e-2 0; 0.106781004451802107743550139118e-1 -0.311642249389153620392187874115e-1 0.150951614734538115101380069557e-1 0.241589449213721869152620314805e-1 -0.249572219856036346729865425447e-1 0.133999434001422053713025372219e0 -0.762334003440718790010098904135e0 0 0.800000000000000000000000000000e0 -0.200000000000000000000000000000e0 0.380952380952380952380952380952e-1 -0.357142857142857142857142857143e-2;];


% construct the lower portion of Q
for i = 1:8
  for j = 1:12
    Q(n-i+1,n-j+1) = -Q(i,j);
  end
end

% construct the interior of Q
for i=9:n-8
  Q(i,i-4:i+4) = [0.1e1 / 0.280e3 -0.4e1 / 0.105e3 0.1e1 / 0.5e1 -0.4e1 / 0.5e1 0 0.4e1 / 0.5e1 -0.1e1 / 0.5e1 0.4e1 / 0.105e3 -0.1e1 / 0.280e3;];
end