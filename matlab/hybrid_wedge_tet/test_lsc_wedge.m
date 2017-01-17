clearvars -global
clear

for ref = 1:5
    switch ref
        case 1
            wedge_convergence1
        case 2
            wedge_convergence2
        case 3
            wedge_convergence3
        case 4
            wedge_convergence4
        case 5
            wedge_convergence5
    end
    
    h = 1.5*(.5.^(ref-1));
    
%     % random perturbations
%     ids = abs(abs(VX)-1) > 1e-8;
%     VX(ids) = VX(ids) + h/4*randn(size(VX(ids)));
%     ids = abs(abs(VY)-1) > 1e-8;
%     VY(ids) = VY(ids) + h/4*randn(size(VY(ids)));
%     ids = abs(abs(VZ)-1) > 1e-8;
%     VZ(ids) = VZ(ids) + h/4*randn(size(VZ(ids)));
    
    % % arnold type mesh
    nlev = 2^(ref-1);
    t = linspace(-1,1,nlev+1);
    for lev = 2:nlev
        ids = find(abs(VZ-t(lev))<1e-8);
        
        [VYlev yids] = sort(VY(ids));
        VXlev = reshape(VX(ids(yids)),nlev+1,nlev+1);
        xids = zeros(nlev+1);
        for j = 1:nlev+1
            [~,xidj] = sort(VXlev(:,j));
            xids(:,j) = xidj + (j-1)*(nlev+1);
        end
        xids = xids(:);       
        p = ids(yids(xids));
        p = reshape(p,nlev+1,nlev+1);
        toggle = [1 0 -1 0]; toggle = toggle(mod(lev-2,4)+1); % switch between add, don't move, subtract
            
        for line = 1:nlev+1
            perturb = repmat([-1 1],1,floor((nlev+1)/2));
            if mod(length(VZ(p)),2)==1
                perturb = [perturb -perturb(end)];
            end
            if mod(line,2)==1                                
                VZ(p(:,line)) = VZ(p(:,line)) + toggle*h*perturb(:);
            else
                VZ(p(:,line)) = VZ(p(:,line)) - toggle*h*perturb(:);
            end
        end
    end
    
    K = size(EToV,1);
    
    filename = sprintf('arnold_warped_%d',ref);
    write_gmsh_file(filename,VX,VY,VZ,EToV);
    
    N = 3;
    
    if 0
        clf
        hold on
        for e = 1:K
            ids = [1 2 3 1 4 5 2 1 3 6 4 1 3 2 5 6];
            plot3(VX(EToV(e,ids)),VY(EToV(e,ids)),VZ(EToV(e,ids)),'linewidth',2)
        end
        keyboard
    end
    
    [r s t w] = wedge_cubature(N);
    
    Nq = length(r);
    x = zeros(Nq,K); y = zeros(Nq,K); z = zeros(Nq,K); J = zeros(Nq,K);
    for e = 1:K
        v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
        [xqe,yqe,zqe,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jqe] = ...
            wedge_geom_factors(VX(v),VY(v),VZ(v),r,s,t);
        x(:,e) = xqe; y(:,e) = yqe; z(:,e) = zqe; J(:,e) = Jqe;
        if (mod(e,round(K/10))==0)
            disp(sprintf('on elem %d out of %d\n',e,K));
        end
    end
    
    Vq = wedge_basis(N,r,s,t);
    M = Vq'*diag(w)*Vq;
    
    a = 1;
    W = (2*a-1)/2*pi;
    f = @(x,y,z) cos(W*x).*cos(W*y).*cos(W*z);
    
    c_lsc = Vq'*diag(w)*(sqrt(J).*f(x,y,z));
    ulsc = (Vq*c_lsc)./sqrt(J);
    
    diff = ulsc - f(x,y,z);
    err2 = diag(w)*(J.*diff.^2);
    errlsc(ref) = sqrt(sum(err2(:)));
    
    for e = 1:K
        M = Vq'*diag(w.*J(:,e))*Vq;
        c(:,e) = M\(Vq'*(w.*J(:,e).*f(x(:,e),y(:,e),z(:,e))));
    end
    diff = Vq*c - f(x,y,z);
    err2 = diag(w)*(J.*diff.^2);
    err(ref) = sqrt(sum(err2(:)));
    
end
%%

h = .5.^(1:length(errlsc));
loglog(h,errlsc,'o-');hold on

start = 3;
h = h(start:end);
errlsc = errlsc(start:end);
C = [h(:).^0 log(h(:))]\log(errlsc(:));
C(2)

h = .5.^(1:length(err));
loglog(h,err,'o-');hold on

h = h(start:end);
err = err(start:end);
C = [h(:).^0 log(h(:))]\log(err(:));
C(2)






