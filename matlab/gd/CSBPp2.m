function [H,Q,X]=CSBPp2(n)
%n = number of nodes in the nodal distribution

Q = zeros(n,n);
H = eye(n,n);
X = zeros(n,1); 

%construct the nodal distribution
dx = 2/(n-1);

for i = 1:n
  X(i) = -1+(i-1)*dx; 
end

Hdiag = [0.17e2 / 0.48e2 0.59e2 / 0.48e2 0.43e2 / 0.48e2 0.49e2 / 0.48e2];

% construct the H matrix_type
for i = 1:4
H(i,i)=Hdiag(i);
H(n-i+1,n-i+1) = Hdiag(i);
end

% add the mesh spacing
H = dx*H;

% construct the upper portion of Q
Q(1:4,1:6) = [-0.1e1 / 0.2e1 0.59e2 / 0.96e2 -0.1e1 / 0.12e2 -0.1e1 / 0.32e2 0 0; -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2 0 0 0; 0.1e1 / 0.12e2 -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2 -0.1e1 / 0.12e2 0; 0.1e1 / 0.32e2 0 -0.59e2 / 0.96e2 0 0.2e1 / 0.3e1 -0.1e1 / 0.12e2;];

% construct the lower portion of Q
for i = 1:4
  for j = 1:6
    Q(n-i+1,n-j+1) = -Q(i,j);
  end
end

% construct the interior of Q
for i=5:n-4
  Q(i,i-2:i+2) = [0.1e1 / 0.12e2 -0.2e1 / 0.3e1 0 0.2e1 / 0.3e1 -0.1e1 / 0.12e2;];
end