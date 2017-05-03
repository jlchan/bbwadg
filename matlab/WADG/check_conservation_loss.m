clear
Globals2D

N = 3
avec = [1e-1 1e-2 1e-3 1e-4];
% avec = 1e-4; % smooth case
for a = avec
    sk = 1;
    for K1D = [4 8 16 32]
        if K1D==1
            [VX VY] = Nodes2D(1); K = 1; EToV = 1:3;
        else
            [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
        end
        
        StartUp2D;
        
        [rq sq w] = Cubature2D(40); % integrate u*v*c
        Vq = Vandermonde2D(N,rq,sq)/V;
        xq = Vq*x; yq = Vq*y;
                
        cfun = @(x,y) 1 + .5*sin(a*pi*x).*sin(a*pi*y); % rough c
%         cfun = @(x,y) 1 + sqrt((x+.1).^2 + (y+.1).^2 + a);
        cfun = @(x,y) exp(x+y);
        
        cq = cfun(xq,yq);
        Mref = Vq'*diag(w)*Vq;
        
        u = @(x,y) exp(x+y);
        u = @(x,y) 1 + sqrt(x.^2 + y.^2 + a) + x;
        %         u = @(x,y) sin(pi*x).*sin(pi*y); % smooth u
        %     a = 16; u = @(x,y) sin(a*pi*x).*sin(a*pi*y); % rough u
%             u = @(x,y) (x > 2*y);%+ (x > -2*y);% + (y>0) + (x > 0);
        %     keyboard
        
        uq = u(xq,yq);
        b = (Vq'*(diag(w)*uq));
        uproj = Mref\b;
        
        %     plot3(x,y,uproj,'o');keyboard
        
        local_cons_err = 0.0;
        fixed_cons_err = 0.0;
        proj_err = 0.0;
        proj_err_weight = 0.0;
        proj_err_cons = 0.0;
        proj_diff = 0.0;
        vmax = 0;
        for e = 1:K
            Minvcu = J(1,e)*Vq'*(w.*1./cq(:,e).*uq(:,e)); % true 1/c mass matrix
            MinvMcMu = J(1,e)*Mref*((Vq'*diag(w.*cq(:,e))*Vq)\(Mref*uproj(:,e))); % weighted DG mass matrix
            
            Minvc = J(1,e)*Vq'*diag(w.*1./cq(:,e))*Vq; % true 1/c mass matrix
            MinvMcM = J(1,e)*Mref*((Vq'*diag(w.*cq(:,e))*Vq)\Mref); % weighted DG mass matrix
%             eK = r(:).^0;
            eK = ones(Np,1); % constant in nodal basis
            v = (MinvMcM-Minvc)*eK;
            alpha = 0;
            if sqrt(abs(v'*eK)) > (1e-8)*norm(v)
                alpha = 1/(v'*eK);
            end
            vmax = max(vmax,norm(v));
            Mcons = MinvMcM - alpha*(v*v');
            
            local_cons_err = local_cons_err + abs(eK'*(Minvcu - MinvMcMu));
            fixed_cons_err = fixed_cons_err + abs(eK'*(Mcons*uproj(:,e) - Minvc*uproj(:,e)));
            %         local_cons_err = local_cons_err + norm(eK'*(Minvc - MinvMcM));
            %         fixed_cons_err = fixed_cons_err + norm(eK'*(Mcons - Minvc));
            
            bK = (Vq'*(J(1,e)*w.*uq(:,e)));
            
            % (u/c,v) = (p,v) -> u/c ~ p
            ucproj = Minvc\bK;
            diff = abs((Vq*ucproj)./cq(:,e)-uq(:,e));
            proj_err = proj_err + sum(sum(J(1,e)*w.*diff.^2));
            
            % (cu,v) = (p,v) -> u/c ~ p
            ucproj = MinvMcM\bK;
            diff = abs((Vq*ucproj)./cq(:,e)-uq(:,e));
            proj_err_weight = proj_err_weight + sum(sum(J(1,e)*w.*diff.^2));
            
            % (cu,v) = (p,v) -> u/c ~ p
            ucproj = Mcons\bK;
            diff = abs((Vq*ucproj)./cq(:,e)-uq(:,e));
            proj_err_cons = proj_err_cons + sum(sum(J(1,e)*w.*diff.^2));
            
            % diff
            diff = Vq*(Minvc\bK - MinvMcM\bK);
            proj_diff = proj_diff + sum(sum(J(1,e)*w.*diff.^2));
            
        end
        vnorm(sk) = vmax;
        err_cons(sk) = local_cons_err; %abs(e(:)'*Minvc*u(:) - e(:)'*MinvMcM*u(:));
        err_cons_fixed(sk) = fixed_cons_err; %abs(e(:)'*Minvc*u(:) - e(:)'*Mcons*u(:));
        
        err_proj(sk) = sqrt(proj_err);
        err_proj_weight(sk) = sqrt(proj_err_weight);
        err_proj_weight_cons(sk) = sqrt(proj_err_cons);
        diff_proj(sk) = sqrt(proj_diff);
        sk = sk + 1;
    end
    
    h = .5.^(1:length(err_cons));
    h = h(:);
    err_cons = err_cons(:);
    err_proj = err_proj(:);
    err_proj_weight = err_proj_weight(:);
    err_proj_weight_cons = err_proj_weight_cons(:);
    diff_proj = diff_proj(:);
    
    figure(1)
    loglog(h,err_cons,'o-')
    hold on;
    loglog(h,abs(err_cons_fixed)+ eps,'s-')
    % loglog(h,vnorm+eps,'s-')
    legend('Conservation error','Fixed conservation error')
    
    ids = 2:length(h);
    fit = [log(h(ids)) ones(size(h(ids)))]\log(err_cons(ids));
    title(sprintf('Order %d: rate = %f\n',N,fit(1)))
    cons_rate = fit(1);
    
    figure(2)
    loglog(h,err_proj,'o-')
    hold on;
    loglog(h,err_proj_weight,'s-')
    loglog(h,err_proj_weight_cons,'x-')
    % loglog(h,diff_proj,'d-')
    
    fit = [log(h(ids)) ones(size(h(ids)))]\log(err_proj(ids));
    fitdiff = [log(h(ids)) ones(size(h(ids)))]\log(diff_proj(ids));
    title(sprintf('Order %d: rate = %f\n',N,fit(1)))
    
    % hold on;loglog(h,vnorm + eps,'s-')
    legend('Projection error','Weighted projection error','Conserved projection error')
    %%
    disp(sprintf('N = %d, cons error rate = %f\n',N,cons_rate))
    format shorte
    % fprintf('\\hline\n')
    disp('L2 errors')
    fprintf('%4.4e & ',err_proj(1:end-1)); fprintf('%4.4e &%f \\\\',err_proj(end),compute_rate(err_proj)); fprintf('\n')
    ids = 3:length(err_proj_weight);
    fprintf('%4.4e & ',err_proj_weight(1:end-1)); fprintf('%4.4e &%f \\\\',err_proj_weight(end),compute_rate(err_proj_weight,ids)); fprintf('\n')
    fprintf('%4.4e & ',err_proj_weight_cons(1:end-1)); fprintf('%4.4e &%f \\\\',err_proj_weight_cons(end),compute_rate(err_proj_weight_cons)); fprintf('\n')
    
    % fprintf('%4.4e & ',err_proj(1:end-1)); fprintf('%4.4e  \\\\',err_proj(end)); fprintf('\n')
    % fprintf('%4.4e & ',err_proj_weight(1:end-1)); fprintf('%4.4e  \\\\',err_proj_weight(end)); fprintf('\n')
    % fprintf('%4.4e & ',err_proj_weight_cons(1:end-1)); fprintf('%4.4e  \\\\',err_proj_weight_cons(end)); fprintf('\n')
    disp('conserv errors')
    fprintf('%4.4e & ',err_cons(1:end-1)); fprintf('%4.4e &%f \\\\',err_cons(end),cons_rate); fprintf('\n')
    fprintf('%4.4e & ',err_cons_fixed(1:end-1)); fprintf('%4.4e & \\\\',err_cons_fixed(end)); fprintf('\n')
    % fprintf('\\hline\n')
    format
    
    r1 = compute_rate(err_proj,2:length(err_proj));
    r2 = compute_rate(err_proj_weight,2:length(err_proj));
    r3 = compute_rate(err_proj_weight_cons,2:length(err_proj));
    [r1 r2 r3]
    
    % print_pgf_coordinates(h,diff_proj)
    disp('L^2 errors')
    print_pgf_coordinates(h,err_proj)
    print_pgf_coordinates(h,err_proj_weight)
    print_pgf_coordinates(h,err_proj_weight_cons)
    disp('conservation errors')
    print_pgf_coordinates(h,err_cons)
    print_pgf_coordinates(h,err_cons_fixed);
    
end

