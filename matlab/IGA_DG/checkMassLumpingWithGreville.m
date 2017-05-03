clear

NB = 3;
for Ksub = ceil(NB/2):4*NB+4
    
    N = NB+Ksub-1;
%     r = JacobiGQ(0,0,N);
    [Nv, VX, Ksub, EToV] = MeshGen1D(-1,1,Ksub);
    
    % squeeze control pts in to reduce CFL restriction    
    re = linspace(-1,1,NB+Ksub)';
    reKsub = linspace(-1,1,Ksub+1)';
    Ve = bspline_basismatrix(NB+1, [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)], reKsub);
    VX = Ve*re; VX = VX(:)';
    
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r(i) = mean(t((i+1):(i+NB))); % greville
    end    
%     r = JacobiGL(0,0,N);
%     plot(t,t*0,'*');hold on
%     plot(rGL,rGL*0,'o');
%     plot(r,r*0,'x');return
        
    [BVDM M Dr] = bsplineVDM(NB,Ksub,r); % VDM for interp, mass, M\S
    
    % imagesc(M*Dr)
    Sr = M*Dr;
    Sr(abs(Sr)<1e-8) = 0;
    
    M(abs(M)<1e-8) =0;
    
    invV = inv(BVDM);
    M = ((M/BVDM)'/BVDM)';
    plot(Ksub,cond(diag(1./sum(M,2))*M),'x');hold on % Greville = h-indep?
end


