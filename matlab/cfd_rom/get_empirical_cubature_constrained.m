% C,d = constraint matrices

function [w id] = get_empirical_cubature_constrained(J,wq,C,d,tol,maxpts)

b = J'*wq; % "exact" integrals
id = [];
r = b;
bnorm = norm(b);
list = 1:size(J,1);
while norm(r)/bnorm > tol && length(id) < maxpts && length(id) < size(J,1)
    
    % greedy pt search
    %[~,idi] = max((J*r)./(norm(r).*sqrt(sum(J.^2,2))));
    [~,idi] = max((J(list,:)*r)./(norm(r).*sqrt(sum(J(list,:).^2,2))));
    id = [id; list(idi)];
    
    % solve for weight, update res
    %w = J(id,:)'\b;
    A = [2*J(id,:)*J(id,:)' C(:,id)';
        C(:,id) zeros(size(C,1))];
    f = [2*J(id,:)*b; d];    
    W = pinv(A)*f;
    w = W(1:size(id)); 

    r = b-J(id,:)'*w;    
    
    keyboard
    
    % remove row
    list(idi) = [];
    
    if mod(length(id),10)==0
        fprintf('%d empirical cubature pts, residual = %g\n',length(id),norm(r)/bnorm)
    end
end

end