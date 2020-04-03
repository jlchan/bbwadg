function [w id] = get_empirical_cubature(J,wq,tol,maxpts)

% options = optimoptions('lsqlin');
% options.FunctionTolerance = 1e-7;

b = J'*wq; % "exact" integrals
id = [];
r = b; % init residual
w = [];
bnorm = norm(b);

list = 1:size(J,1);
while norm(r)/bnorm > tol & length(id) < maxpts
    
    % greedy pt search
    [~,idi] = max((J(list,:)*r)./(norm(r).*sqrt(sum(J(list,:).^2,2))));    
    id = [id; list(idi)];
    
    % solve for weight, update res
    w = J(id,:)'\b;
    r = b-J(id,:)'*w;
    
    % if weights negative, recompute with constraints 
    if min(w(:))<0 
        disp('running nonlinear lsq opt');
        
        Aeq = []; beq = [];
        w = lsqlin(J(id,:)',b,[],[],Aeq,beq,zeros(size(w)),sum(wq)*ones(size(w)),w); % upper bound = total vol
        r = b-J(id,:)'*w;
    end
    
    % remove index
    list(idi) = [];
    
    if mod(length(id),10)==0
        fprintf('%d empirical cubature pts, residual = %g\n',length(id),norm(r)/bnorm)
    end
end

fprintf('%d empirical cubature nodes, residual for emp cubature = %g\n',length(w),norm(J(id,:)'*w - J'*wq)/norm(J'*wq))

% if min(w(:))<0
%     disp('running nonlinear lsq opt');
%     w = lsqlin(J(id,:)',b,[],[],[],[],zeros(size(w)),sum(wq)*ones(size(w)),w); % upper bound = total vol
%     r = b-J(id,:)'*w;
%
%     fprintf('After positivity correction: %d empirical cubature nodes, residual for emp cubature = %g\n',length(w),norm(J(id,:)'*w - J'*wq)/norm(J'*wq))
% end





