function id = get_empirical_cubature(Vtarget,b,tol,maxpts)

J = Vtarget;
id = [];
r = b;
list = 1:size(J,1);
while norm(r)/norm(b) > tol && length(id) < maxpts
    
    % greedy pt search
    [~,idi] = max((J*r)./(norm(r).*sqrt(sum(J.^2,2))));
    id = [id; list(idi)];
    
    % solve for weight, update res
    w = Vtarget(id,:)'\b;
    r = b-Vtarget(id,:)'*w;
%     r = Vtarget(id,:)'*w-b; % wrong!
    
    % remove row
    J(idi,:) = [];
    list(idi) = [];
    
    if mod(length(id),10)==0
        fprintf('%d empirical cubature pts, residual = %g\n',length(id),norm(r)/norm(b))
    end
end

% fprintf('Empirical cubature residual = %g\n',norm(r)/norm(b))
end