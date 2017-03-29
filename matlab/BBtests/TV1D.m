function TV = TV1D(u)

Globals1D

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end

tol = 1e-2; %max(abs(u - repmat(mean(u,1),N+1,1)))
for i = 2:N
    flag1 = abs(u(i,:)-u(i-1,:)) > tol;
    flag2 = abs(u(i+1,:)-u(i,:)) > tol;
%     TV = TV + abs(sign(u(i,:) - u(i-1,:)).*flag1 - sign(u(i+1,:)- u(i,:)).*flag2); % sign variations at each node    
%     TV = TV + abs(u(i+1,:) - 2*u(i,:) + u(i-1,:)); % sign variations at each node    
end
