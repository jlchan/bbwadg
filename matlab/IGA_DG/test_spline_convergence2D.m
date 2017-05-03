clear
Globals2D;

K1D = 1;
useQuadrature = 1;

x0 = .1; y0 = .1;
pex = @(x,y,t) exp(-5^2*((x-x0).^2 + (y-y0).^2));

k = 5;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2);
% pex = @(x,y,t) sin(k*pi*x).*sin(k*pi*y); 
% pex = @(x,y,t) exp(x.*y).*sin(k*pi*x).*sin(k*pi*y); 
% pex = @(x,y,t) exp(x.*y);

a = .125;
% a = 1/64;
% a = 1/512;
        
% figure
for ii = 1:4    
    if ii==1
        K1D = 1;
        Nmax = 16;
        fKsub = @(NB) 1;
        smoothKnots = 0;
    elseif ii==3
        K1D = 1;
        Nmax = 5;
        fKsub = @(NB) NB;
        smoothKnots = 'opt';
    elseif ii==4
        K1D = 1;
        Nmax = 8;
        fKsub = @(NB) NB;
        smoothKnots = 25;
    elseif ii==2
        K1D = 1;
        Nmax = 8;
        fKsub = @(NB) NB;
        smoothKnots = 0;        
    end
    
    clear L2err ndofs
    sk = 1;
    for NB = 2:Nmax
        
        Ksub = fKsub(NB);
        N = NB+Ksub-1;
        dofs = (N+1)^2*K1D^2;
        
        % Read in Mesh
        [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);
        StartUpQuad2D;
        
        % non-affine mappings
        x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
        y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);
        %     x = x + a*(1-y).*(1+y).*(1-x).*(1+x);
        %     y = y + a*(1-x).*(1+x).*(1-y).*(1+y);
        if 0 && NB>7
            x = reshape(x,N+1,N+1);
            y = reshape(y,N+1,N+1);
            rp1D = linspace(-1,1,100)';
            Vp = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
            hold on
            for i = 1:N+1
                plot(Vp*x(:,i),Vp*y(:,i),'k--');                
            end
            for i = 1:N+1
                plot(Vp*x(i,:)',Vp*y(i,:)','k--');                
            end
            plot(x(:),y(:),'ko','markersize',8,'MarkerFaceColor',[.49 1 .63]);hold on            
            
            axis off
            axis tight
%             print(gcf,'-dpng','../talks/UT2017/figs/mapped.png')
            return
        end
        
        % isoparametric - interpolate geofacs
        [rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
        
%         % initialize TP operators
%         rp1D = linspace(-1,1,150);
%         [rp sp] = meshgrid(rp1D);
%         rp = rp(:); sp = sp(:);
%         Vp = Vandermonde2DQuad(N,rp,sp)/V;
%         xp = Vp*x;
%         yp = Vp*y;
        
        [rq1D wq1D] = JacobiGQ(0,0,2*N+2);
        if (Ksub>1)                       
            [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,rq1D,smoothKnots);    
            [rgq wgq] = JacobiGQ(0,0,2*NB+2);            

            h = @(r) repmat(diff(VX(:)'),length(r),1);
            map = @(r) reshape(h(r).*(repmat(r,1,Ksub)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);
            
            rq1D = map(rgq);            
            wq1D = repmat(wgq,1,Ksub).*h(rgq)/2;
            rq1D = rq1D(:); wq1D = wq1D(:);
        end
        
        [r1D w1D] = JacobiGL(0,0,N);
        if Ksub==1
%            rq1D = r1D;  wq1D = w1D;
        end
%         if Ksub > 1            
%             t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
%             for i = 1:N+1
%                 r1D(i) = mean(t((i+1):(i+NB))); % greville
%             end
%         end
        
        [rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
        [wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
        V1D = Vandermonde1D(N,r1D);
        Vq1D = Vandermonde1D(N,rq1D)/V1D;
        Nq = length(wq1D);
        xq = zeros(Nq^2,K);
        yq = zeros(Nq^2,K);
        Jq = zeros(Nq^2,K);
        for e = 1:K
            xq(:,e) = reshape((Vq1D*reshape(x(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
            yq(:,e) = reshape((Vq1D*reshape(y(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
            Jq(:,e) = reshape((Vq1D*reshape(J(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
        end
        wJq = spdiag(wq)*Jq;
        
        if Ksub==1
            V1D = Vandermonde1D(N,r1D);
%             Vp1D = Vandermonde1D(N,rp1D)/V1D;
            Vq1D = Vandermonde1D(N,rq1D)/V1D;
            M1D = Vq1D'*diag(wq1D)*Vq1D;
            invM1D = V1D*V1D';
        else
            [BVDM M1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
%             Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
            Vq1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
            invM1D = inv(M1D);
        end
        Pq1D = invM1D*Vq1D'*diag(wq1D);
        
%         Vq = kron(Vq1D,Vq1D);
        
        err2 = 0;
        err = zeros(Nq^2,K);
        for e = 1:K
            
            p = (Pq1D*reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq))*Pq1D';
%             p = (Vq'*diag(wq.*Jq(:,e))*Vq)\(Vq'* (wq.*Jq(:,e).*pex(xq(:,e),yq(:,e),0))); p = reshape(p,N+1,N+1);            
            pterr = (Vq1D*p)*Vq1D' - reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq);
            err(:,e) = pterr(:);
            errK = reshape(wJq(:,e),Nq,Nq).*(pterr).^2;
            err2 = err2 + sum(errK(:));
        end
        
%         if ii==1 && N==12
%             keyboard
%         end
        L2err(sk) = sqrt(err2);
        ndofs(sk) = dofs;
        sk = sk + 1;
    end
    
    if ii==1
        semilogy(ndofs.^(.5),L2err,'.-','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = %d',Ksub,smoothKnots))
    elseif ii==2        
        semilogy(ndofs.^(.5),L2err,'o--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = opt',Ksub))
    elseif ii==3
        semilogy(ndofs.^(.5),L2err,'x--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = %d',Ksub,smoothKnots))
    else
        semilogy(ndofs.^(.5),L2err,'d--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = %d',Ksub,smoothKnots))
    end
    hold on
    
    if strcmp(smoothKnots,'opt')
%         disp('optknots =====')
        print_pgf_coordinates(ndofs,L2err)
%         disp('=====')
    else
%         disp(sprintf('smooth = %d ======== \n',smoothKnots))
        print_pgf_coordinates(ndofs,L2err)
%         disp('=====')
    end
end
title(sprintf('K = %d',K))
grid on
legend show


%% href 

clear
Globals2D;

K1D = 1;
useQuadrature = 1;

x0 = .1; y0 = .1;
pex = @(x,y,t) exp(-5^2*((x-x0).^2 + (y-y0).^2));

k = 1;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2);
% pex = @(x,y,t) sin(k*pi*x).*sin(k*pi*y); 
% pex = @(x,y,t) exp(x.*y).*sin(k*pi*x).*sin(k*pi*y); 
% pex = @(x,y,t) exp(x.*y);

a = .125;
        
% figure
for ii = 1:3
    ii
    if ii==1
        K1D = 2;
        NB = 4;
        Kvec = [4 8 16 32];
        smoothKnots = 0;            
    elseif ii==2
        K1D = 2;
        NB = 4; 
        Kvec = [4 8 16 32];
        smoothKnots = 'opt';
    elseif ii==3
        K1D = 2;
        NB = 4; 
        Kvec = [4 8 16 32];
        smoothKnots = 50;
    end
    
    clear L2err ndofs
    if 1
        sk = 1;
        for Ksub = Kvec
            
            N = NB+Ksub-1;
            dofs = (N+1)^2*K1D^2;
            
            % Read in Mesh
            [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);
            StartUpQuad2D;
            
            % non-affine mappings
            if a > 1e-10
                x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
                y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);
                %     x = x + a*(1-y).*(1+y).*(1-x).*(1+x);
                %     y = y + a*(1-x).*(1+x).*(1-y).*(1+y);
                %         if NB>6
                %             plot(x,y,'o');return
                %         end
            end
            
            % isoparametric - interpolate geofacs
            [rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
            
            [rq1D wq1D] = JacobiGQ(0,0,N+1);
            if (Ksub>1)
                [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
                [rgq wgq] = JacobiGQ(0,0,NB+1);
                
                h = @(r) repmat(diff(VX(:)'),length(r),1);
                map = @(r) reshape(h(r).*(repmat(r,1,Ksub)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);
                
                rq1D = map(rgq);
                wq1D = repmat(wgq,1,Ksub).*h(rgq)/2;
                rq1D = rq1D(:);
                wq1D = wq1D(:);
            end
            
            [r1D] = JacobiGL(0,0,N);
            
            [rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
            [wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
            V1D = Vandermonde1D(N,r1D);
            Vq1D = Vandermonde1D(N,rq1D)/V1D;
            Nq = length(wq1D);
            xq = zeros(Nq^2,K);
            yq = zeros(Nq^2,K);
            Jq = zeros(Nq^2,K);
            for e = 1:K
                xq(:,e) = reshape((Vq1D*reshape(x(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
                yq(:,e) = reshape((Vq1D*reshape(y(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
                Jq(:,e) = reshape((Vq1D*reshape(J(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
            end
            wJq = spdiag(wq)*Jq;
            
            
            if Ksub==1
                V1D = Vandermonde1D(N,r1D);
                D1D = GradVandermonde1D(N,r1D)/V1D;
                %             Vp1D = Vandermonde1D(N,rp1D)/V1D;
                %             Vq1D = Vandermonde1D(N,rq1D)/V1D;
                M1D = Vq1D'*diag(wq1D)*Vq1D;
                invM1D = V1D*V1D';
            else
                [BVDM M1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
                %             Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
                Vq1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
                invM1D = inv(M1D);
            end
            Pq1D = invM1D*(Vq1D'*diag(wq1D));
            
            %         Vq = kron(Vq1D,Vq1D);
            
            err2 = 0;
            err = zeros(Nq^2,K);
            for e = 1:K
                
                p = (Pq1D*reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq))*Pq1D';
                %             p = (Vq'*diag(wq.*Jq(:,e))*Vq)\(Vq'* (wq.*Jq(:,e).*pex(xq(:,e),yq(:,e),0))); p = reshape(p,N+1,N+1);
                pterr = (Vq1D*p)*Vq1D' - reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq);
                err(:,e) = pterr(:);
                errK = reshape(wJq(:,e),Nq,Nq).*(pterr).^2;
                err2 = err2 + sum(errK(:));
            end
            
            %         if ii==1 && N==12
            %             keyboard
            %         end
            L2err(sk) = sqrt(err2);
            ndofs(sk) = dofs;
            sk = sk + 1;
        end
    end
    
    if ii==1
        semilogy(ndofs.^(.5),L2err,'o--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = opt',Ksub))
    elseif ii==2
%         keyboard
        semilogy(ndofs.^(.5),L2err,'x--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = %d',Ksub,smoothKnots))
    else
        semilogy(ndofs.^(.5),L2err,'d--','linewidth',2,'DisplayName',sprintf('Ksub = %d, smoothKnots = %d',Ksub,smoothKnots))
    end
    hold on
    
    if ii==1
%         disp('optknots =====')
        print_pgf_coordinates(ndofs,L2err)
%         disp('=====')
    else
%         disp(sprintf('smooth = %d ======== \n',smoothKnots))
        print_pgf_coordinates(ndofs,L2err)
%         disp('=====')
    end
end
title(sprintf('K = %d',K))
grid on
legend show
