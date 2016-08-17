clear

for useBern = 1
    
    Globals1D
    
    for K = 3; %2.^[4:5]
        
        sk = 1;
        for N = 5 %1:4
            
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
                        
%             KK = Dr'*MK*Dr;
%             KK(1,:) = 0; KK(:,1) = 0; KK(1,1) = 1;
%             KK(end,:) = 0; KK(:,end) = 0; KK(end,end) = 1;            
%             keyboard
            
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
                        
            b = sum(M,2);
            
            % BCs
            Ks(1,:) = 0; Ks(:,1) = 0; Ks(1,1) = 1;
            Ks(end,:) = 0; Ks(:,end) = 0; Ks(end,end) = 1;
            b(1) = 0; b(end) = 0;
            
            % preconditioner matrix                       
            NpK = size(Ks,2); % global dofs
            Kpre = 2*diag(ones(NpK,1)) - diag(ones(NpK-1,1),1) - diag(ones(NpK-1,1),-1);
            Kpre = (NpK-1)*Kpre; 
            Kpre(1,:) = 0; Kpre(:,1) = 0; Kpre(1,1) = 1;
            Kpre(end,:) = 0; Kpre(:,end) = 0; Kpre(end,end) = 1;            
            
            iK = R*inv(Ks);
            for i = N+1:2*N+1%max(gid(:))
                e = zeros(max(gid(:)),1); e(i) = 1;                
                g = reshape(iK*e,N+1,K);
                vv = Vp*g; vv = vv/max(abs(vv(:)));
%                 plot(xp,vv);hold on
                
                vv = vv(:,2);
                v0 = vv(1)*(1-rp)/2 + vv(end)*(1+rp)/2;                
                vv = vv-v0; 
                if (i>N+1) && (i<2*N+1)
                    vv = vv/max(abs(vv(:)));
                end
                plot(rp,vv);hold on
%                                 
%                 rSEM = JacobiGL(0,0,N);
%                 hold on;plot(rSEM,rSEM*0,'o','markersize',12)
            end            
            
            keyboard
                        
            kappa(sk) = cond(Kpre\Ks);
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

