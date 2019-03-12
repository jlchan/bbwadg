function id = get_empirical_cubature(J,b,tol,maxpts)

id = [];
r = b;
bnorm = norm(b);
list = 1:size(J,1);
while norm(r)/bnorm > tol && length(id) < maxpts
    
    % greedy pt search
    %[~,idi] = max((J*r)./(norm(r).*sqrt(sum(J.^2,2))));
    [~,idi] = max((J(list,:)*r)./(norm(r).*sqrt(sum(J(list,:).^2,2))));
    id = [id; list(idi)];
    
    % solve for weight, update res
    w = J(id,:)'\b;
    r = b-J(id,:)'*w;    
    
    % remove row
    list(idi) = [];
    
    if mod(length(id),10)==0
        fprintf('%d empirical cubature pts, residual = %g\n',length(id),norm(r)/bnorm)
    end
end

end