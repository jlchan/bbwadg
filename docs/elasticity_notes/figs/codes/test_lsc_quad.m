% Driver script for solving the 2D Poisson equation
clear
useQuads=1; mypath
Globals2D

N = 3;

% random pertrubations
ax = 0; ay = ax;  % for uniform
% ax = 1.5/N; ay = ax;  % for perturbed uniform
% ax = 1/N; ay = 1/N; % for skewing Arnold
delta_skew = 0;

deltaCurv = .25;
plot_meshes = 2;
Nlev = 4;

% reference nodes
r1 = [-1 1 1 -1]; s1 = [-1 -1 1 1];
r1D = JacobiGL(0,0,N);
[r s] = meshgrid(r1D);
r = r(:); s = s(:);
Nq = 3*N;
[rq1D w1D] = JacobiGQ(0,0,Nq);
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(w1D); w = wr(:).*ws(:);

V1 = Vandermonde2D(1,r1,s1);
V1N = Vandermonde2D(1,r,s)/V1;
V = Vandermonde2D(N,r,s);
Vq = Vandermonde2D(N,rq,sq);
[Vr, Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/V;  Ds = Vs/V;
[Vrq, Vsq] = GradVandermonde2D(N,rq,sq);
Drq = Vrq/V;  Dsq = Vsq/V;

rp1D = linspace(-1,1,150);
V1D = Vandermonde1D(N,r1D);
Vp1D = Vandermonde1D(N,rp1D)/V1D;

for lev=1:Nlev
    lev
    % build arnold mesh
    Nx = 2^lev+1;
    x1d = linspace(-1,1,Nx);
    [VYQ,VXQ] = meshgrid(x1d);
    
    EToVQ = [];
    k = 1;
    for m=1:Nx-1
        for n=1:Nx-1
            EToVQ(k,:) = [n+(m-1)*Nx, n+1+(m-1)*Nx, n+1+m*Nx, n+m*Nx];
            k = k+1;
        end
    end
    hlev = 1/Nx;
    delta = delta_skew*hlev;
    VYQ(2:2:end,2:2:end) = VYQ(2:2:end,2:2:end) + delta;
    VYQ(1:2:end,2:2:end) = VYQ(1:2:end,2:2:end) - delta;
    
    EToV = EToVQ; K = size(EToV,1);
    VX = VXQ(:)'; VY = VYQ(:)';
    
    %         perturb = @(y) .5*sin(Nx*pi*y)/Nx^2;
    %         perturb = @(y) 2*y/Nx^2;
    %         perturb = @(y) 0*y/Nx^2;
    %         ids = abs(abs(VX)-1)>1e-8;
    %         VX(ids) = VX(ids) + perturb(VY(ids));
    %         ids = abs(abs(VY)-1)>1e-8;
    %         VY(ids) = VY(ids) + perturb(VX(ids));
    
    [EToEQ, EToFQ] = tiConnectQuad2D(EToV);
    EToE = EToEQ; EToF = EToFQ;
    StartUpQuad2D;
    
    % unique vertex nodes: search nodes for vertex matches
    Vmask = [1 N+1 (N+1)^2 (N+1)^2-(N+1)];
    for e = 1:K
        for v = 1:4
            idv = find(abs(x - x(Vmask(v),e))<1e-8 & abs(y - y(Vmask(v),e))<1e-8);
            vmapV{e,v} = idv; % list of vertices connected to this vertex
        end
    end
    
    vmapM = reshape(vmapM,Nfp*Nfaces,K);
    vmapP = reshape(vmapP,Nfp*Nfaces,K);
    xb = x(vmapB);
    yb = y(vmapB);
    for e = 1:K
        xe = x(:,e);
        ye = y(:,e);
        
        dx = randn(size(xe)); dx = dx/max(abs(dx(:)));
        x(:,e) = x(:,e) + ax*dx/Nx^2;
        dy = randn(size(ye)); dy = dy/max(abs(dy(:)));
        y(:,e) = y(:,e) + ay*dy/Nx^2;
        
        vmapMe = reshape(vmapM(:,e),Nfp,Nfaces);
        vmapPe = reshape(vmapP(:,e),Nfp,Nfaces);
        
        % carry perturbation over to neighbors
        for f = 1:4 % nfaces = 4
            x(vmapPe(:,f)) = x(vmapMe(:,f));
            y(vmapPe(:,f)) = y(vmapMe(:,f));
        end
        for v = 1:4
            x(vmapV{e,v}) = x(Vmask(v),e);
            y(vmapV{e,v}) = y(Vmask(v),e);
        end
    end
    x(vmapB) = xb;
    y(vmapB) = yb;
    
    % sinusoidal perturbations
    for e = 1:K
        col = mod(e-1,Nx-1) + 1;
        row = floor((e-1)/(Nx-1)) + 1;
        dy = deltaCurv*(-1)^(row+1)*cos((Nx-1)/2*pi*(x(:,e)+1))/Nx;
        if (mod(row,2)==1) % if odd row, use top face. if even, do nothing
            f = 3;
            y(Fmask(:,f),e) = y(Fmask(:,f),e) + dy(Fmask(:,f));
            
            % match
            vmapMe = reshape(vmapM(:,e),Nfp,Nfaces);
            vmapPe = reshape(vmapP(:,e),Nfp,Nfaces);
            y(vmapPe(:,f)) = y(vmapMe(:,f));
            
            % blend displacement
            s1D = JacobiGL(0,0,N);
            for slice = 1:N
                ids = (1:N+1) + (slice-1)*(N+1);
                scale = (1+s1D(slice))/2;
                y(ids,e) = y(ids,e) + dy(Fmask(:,f))*scale;
            end
            
            % blend nbr displacement
            enbr = EToE(e,f);  fnbr = 1;
            for slice = 2:N+1
                ids = (1:N+1) + (slice-1)*(N+1);
                scale = (1 - s1D(slice))/2;
                y(ids,enbr) = y(ids,enbr) + dy(Fmask(:,fnbr))*scale;
            end
        end
    end
    
    % recompute J
    for e = 1:K
        xe = x(:,e);
        ye = y(:,e);
        %[rx,sx,ry,sy,Je] = GeometricFactorsQuad2D(xe,ye,Dr,Ds);
%         J(:,e) = Je;
% Jq(:,e) = (Vq/V)*Je;
% [rx,sx,ry,sy,Je] = GeometricFactorsQuad2D(xe,ye,Drq,Dsq);
% Jq(:,e) = -Je;
        xq(:,e) = (Vq/V)*xe;
        yq(:,e) = (Vq/V)*ye;
        
    end
    
    Jq = zeros(size(xq));
    [r2 s2] = Nodes2D(2*N); [r2 s2] = xytors(r2,s2);
    V2 = Vandermonde2D(N,r2,s2);
    x2 = (V2/V)*x;    y2 = (V2/V)*y;
    V2 = Vandermonde2D(2*N,r2,s2);
    [Vr2 Vs2] = GradVandermonde2D(2*N,r2,s2);
    Dr2 = Vr2/V2;    Ds2 = Vs2/V2;
    Vq2 = Vandermonde2D(2*N,rq,sq);
    [rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x2,y2,Dr2,Ds2); Jq = (Vq2/V2)*J;   % high order polynom
    
    
    
    if min(Jq(:))<1e-8
        plot(x,y,'o')
        %         plot(xq,yq,'.')
        keyboard
    end
    
    if lev == plot_meshes
        figure(lev)
        hold on;
        for e = 1:K
            xe = Vp1D*reshape(x(:,e),N+1,N+1);
            ye = Vp1D*reshape(y(:,e),N+1,N+1);
            plot(xe,ye,'k--','linewidth',1);
            xe = reshape(x(:,e),N+1,N+1)*Vp1D';
            ye = reshape(y(:,e),N+1,N+1)*Vp1D';
            plot(xe',ye','k--','linewidth',1);
            
            for f = 1:4
                xe = Vp1D*x(Fmask(:,f),e);
                ye = Vp1D*y(Fmask(:,f),e);
                plot(xe,ye,'k-','linewidth',2);
            end
        end
        
        plot(x,y,'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
        axis off
        axis equal
        %             axis tight
        return;
    end
    
    % target function
    alpha = 1;
    beta = 1;
    %f = @(x,y) sin(alpha*pi*x).*sin(beta*pi*y);
    f = @(x,y) exp(alpha*x+beta*y);
    %     f = @(x,y) x.^N.*y + x.*y.^N;
    fq = f(xq,yq);
    
    errl2 = 0;
    errlsc = 0;
    errproj = 0;
    Mref = (Vq'*diag(w)*Vq); % refernece M
    for e=1:K
        
        % L2 projection
        u = (Vq'*diag(w.*Jq(:,e))*Vq) \ (Vq'*(w.*Jq(:,e).*fq(:,e)));
        err = sum(sum(w.*Jq(:,e).*(Vq*u - fq(:,e)).^2));
        uproj = u;
        errl2 = errl2 + err;
        
        % LSC-DG
        u = Mref \ (Vq'*(w.*sqrt(Jq(:,e)).*fq(:,e)));
        err = sum(sum(w.*Jq(:,e).*((Vq*u)./sqrt(Jq(:,e)) - fq(:,e)).^2));
        errlsc = errlsc + err;
        
        % weighted projection
        MinvJ = (Vq'*diag(w./Jq(:,e))*Vq);
        Mw = Mref*(MinvJ\Mref);
        u = Mw \ (Vq'*(w.*Jq(:,e).*fq(:,e)));
        err = sum(sum(w.*Jq(:,e).*(Vq*u - fq(:,e)).^2));
        %         err = sum(sum(w.*Jq(:,e).*(Vq*u - Vq*uproj).^2)); % consistency error
        errproj = errproj + err;
        
    end
    errl2save(lev) = sqrt(errl2);
    errlscsave(lev) = sqrt(errlsc);
    errprojsave(lev) = sqrt(errproj);
    minJ(lev) = min(J(:));
    %     err(lev) =  max(max(abs(ulsc-quQ)))
    h(lev) = 2/Nx;
end

h = .5.^(1:lev);

if plot_meshes>0
    %print(figure(1),'-dpng','../arnoldmesh1.png')
    %print(figure(2),'-dpng','../arnoldmesh2.png')
    %     print(figure(1),'-dpng','../randunif1.png')
    %     print(figure(2),'-dpng','../randunif2.png')
    disp('printing meshes to file')
    %     print(gcf,'-dpng','../curvarnold1.png')
    
    %print(gcf,'-dpng','../curvarnold3.png')
    return
end


print_pgf_coordinates(h,errl2save)
print_pgf_coordinates(h,errlscsave)
print_pgf_coordinates(h,errprojsave)


figure
loglog(h,errl2save, 'r-*');
hold on
loglog(h,errlscsave, 'b-o');
loglog(h,errprojsave, 'ks-');
% loglog(h,minJ, 'm^-');

len = 1;
h = h(end-len:end);
errl2save = errl2save(end-len:end);
errprojsave = errprojsave(end-len:end);
minJ = minJ(end-len:end);
fit1 = [log(h(:)) ones(size(h(:)))]\log(errl2save(:));
fit2 = [log(h(:)) ones(size(h(:)))]\log(errprojsave(:));
fitJ = [log(h(:)) ones(size(h(:)))]\log(minJ(:));

% title(sprintf('L2 proj rate = %f vs projected mass = %f. min J = rate %f\n',fit1(1),fit2(1),fitJ(1)))
title(sprintf('L2 proj rate = %f vs projected mass = %f\n',fit1(1),fit2(1)))

legend('L2 projection','LSC-DG','Weighted projection')
set(gca,'fontsize',14)
axis equal

useQuads=0; mypath