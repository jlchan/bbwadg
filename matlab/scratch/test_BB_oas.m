clear

for useBern = 0:1
    
    Globals1D
    
    for K = 4%2.^[4:5]
        
        sk = 1;
        for N = 4%1:25
            
            % Generate simple mesh
            [Nv, VX, K, EToV] = MeshGen1D(0,1,K);
            
            %N,K,size(EToV),size(VX)
            
            StartUp1D
            map = @(r) ones(length(r),1)*VX(va) + .5*(r(:)+1)*(VX(vb)-VX(va));
            h = 1/K;
            
            % deriv ops
            r = JacobiGL(0,0,N);
            [rq w] = JacobiGQ(0,0,N);
            xq = map(rq);
            
            % plotting
            Nplot = 250;
            rp = linspace(-1,1,Nplot); rp = rp(:);
            xp = map(rp);
            
            if useBern
                [V Vr] = bern_basis_1D(N,r); Dr = V\Vr;
                Vq = bern_basis_1D(N,rq);
                Vp = bern_basis_1D(N,rp);
                EK = bern_basis_1D(N,r)\bern_basis_1D(1,r);
            else
                % lumping option for SEM
                rq = JacobiGL(0,0,N); VSEM = Vandermonde1D(N,rq); w = sum(inv(VSEM*VSEM'),2);
                
                V = Vandermonde1D(N,r); Vr = GradVandermonde1D(N,r);  Dr = Vr/V;
                Vq = Vandermonde1D(N,rq)/V;
                Vp = Vandermonde1D(N,rp)/V;
                
                EK = (Vandermonde1D(N,r)/V)\bern_basis_1D(1,r);
                
            end
            
            % mass and stiffness
            MK = Vq'*diag(w)*Vq;
            
            % assembly op
            gid = zeros((N+1),K);
            off = 0;
            for e = 1:K
                gid(:,e) = off + (1:N+1);
                off = off + N;
            end
            R = zeros((N+1)*K,max(gid(:)));
            for i = 1:length(gid(:))
                R(i,gid(i)) = 1;
            end
            
            % assembled global matrices
            Ks = R'*(2/h)*kron(eye(K),Dr'*MK*Dr)*R;
            M = R'*(h/2)*kron(eye(K),MK)*R;
            
            %             Ks = M + Ks; % helmholtz
            b = sum(M,2);
            
            % BCs
            Ks(1,:) = 0; Ks(:,1) = 0; Ks(1,1) = 1;
            Ks(end,:) = 0; Ks(:,end) = 0; Ks(end,end) = 1;
            b(1) = 0; b(end) = 0;
            
            % preconditioner matrix
            Kpre = zeros(size(Ks));
            ids = -1:N+1; % precon nodes
            % assemble left, middle, right precon from assembled mat
            KL = Ks(1:N+2,1:N+2);
            KM = Ks(N:2*N+2,N:2*N+2);
            Kpre(1:N+2,1:N+2) = inv(KL);
            off = N+1;
            for e = 2:K-1
                Kpre(off+ids,off+ids) = Kpre(off+ids,off+ids) + inv(KM);
                off = off + N; % move fwd one element
            end
            Kpre(end-(N+1):end,end-(N+1):end) = Kpre(end-(N+1):end,end-(N+1):end)+inv(KL(N+2:-1:1,N+2:-1:1)); % permute by symmetry
            
            % coarse grid matrix
            Kcoarse = 2*diag(ones(K-1,1)) - diag(ones(K-2,1),1) - diag(ones(K-2,1),-1);
            Kcoarse = (1/h)*Kcoarse;
            Kcoarse = [1 zeros(1,K); zeros(K-1,1) Kcoarse zeros(K-1,1); zeros(1,K) 1];
            
            % coarse grid assembly
            R0 = zeros(2*K,length(VX));
            off = 0;
            for e = 0:K-1
                R0(2*e + (1:2),off + (1:2)) = eye(2);
                off = off + 1;
            end
            
            E = R0'*kron(eye(K),EK')*R*diag(1./sum(R,1)); % transfer from C0 degree N to linear
            % transfer op from degree N to linear with BCs built in
            Ebc = diag(1./sum(E,2))*E;
            Ebc(1,:) = 0; Ebc(1,1) = 1;
            Ebc(end,:) = 0; Ebc(end,end) = 1;
            Kcoarse = E'*inv(Kcoarse)*Ebc;
            
            % a = 1 => add coarse grid
            a = 1; 
            if 0 % plot greens function           
                e = zeros(size(Kpre,2),1);
                e(gid(abs(.5-x(:))<1e-8)) = 1;                
                u = reshape(R*((0*Kcoarse+Kpre)*e),N+1,K); u = u/max(abs(u(:)));
                if useBern
                    plot(xp,Vp*u,'b.-','DisplayName',sprintf('h = %g, N = %d, Bernstein = %d',h,N,useBern));hold on
                else
                    plot(xp,Vp*u,'r.-','DisplayName',sprintf('h = %g, N = %d, Bernstein = %d',h,N,useBern));hold on
                end
                u = reshape(R*(Ks\e),N+1,K); u = u/max(abs(u(:)));
                plot(xp,Vp*u,'--','linewidth',2,'DisplayName',sprintf('h = %g, N = %d, Bernstein = %d',h,N,useBern));hold on
%                 keyboard
            end
            kappa(sk) = cond((a*Kcoarse+Kpre)*Ks);
            sk = sk + 1;
        end
        semilogy(kappa,'o-','DisplayName',sprintf('h = %g, Bernstein = %d',h,useBern));hold on;
        xlabel('Degree N','fontsize',14)
    end
    if useBern
        kappaBern = kappa;
    else
        kappaSEM = kappa;
    end
    
end
legend show

return

u = reshape(R*(Ks\b),N+1,K);

% rg = equispaced pts
rg = linspace(-1,1,N+1); rg = rg(:);
xg = map(rg);

vv = Vp*u;
% plot(rp,Vp);return
hold on
plot(xp,vv,'o-','linewidth',2)
plot(xp,.5*(xp).*(1-xp),'--')

return

