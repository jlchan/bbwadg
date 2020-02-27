Globals2D;

N = 2;

filename = 'Grid/Other/block2.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% K1D = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% This builds the nodal DG stuff
StartUp2D;

Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); 
[rq1D wq1D] = JacobiGQ(0,0,N);
useSBP = 0;

[rq sq wq] = Kubatko2D(N,'Lobatto'); 
[rq1D wq1D] = JacobiGL(0,0,N+1);
useSBP = 1;

Vq = Vandermonde2D(N,rq,sq)/V;

xq = Vq*x;
yq = Vq*y;

rf = [rq1D; -rq1D; -ones(size(rq1D))];
sf = [-ones(size(rq1D)); rq1D; -rq1D];
wf = [wq1D; wq1D; wq1D];
Vf = Vandermonde2D(N,rf,sf)/V;
xf = Vf*x;
yf = Vf*y;



PlotMesh2D
hold on
% for e = 1:K
%     ids = EToV(e,:);
%     text(mean(VX(ids)),mean(VY(ids)),num2str(e));
% end
e = [223 222]; % 209 226]; % center = 223
    
if ~useSBP
    for ee1 = 1
        for ee2 = 2
            for i = 1:length(rq)
                for j = 6%1:length(rq)
                    xx = [xq(i,e(ee1)) xq(j,e(ee2))];
                    yy = [yq(i,e(ee1)) yq(j,e(ee2))];
                    plot(xx,yy,'m:','linewidth',2)
                end
                for j = 1:length(rq)
                    xx = [xq(i,e(ee1)) xq(j,e(ee2))];
                    yy = [yq(i,e(ee1)) yq(j,e(ee2))];
                    plot(xx,yy,'m:')
                end
            end
        end
    end
end
plot(xq(:,e),yq(:,e),'bo','markersize',24,'linewidth',3,'MarkerFaceColor',[.49 1 .63])
fids = (1:length(rq1D))+length(rq1D);
plot(xf(fids,e(1)),yf(fids,e(1)),'rs','markersize',24,'linewidth',3,'MarkerFaceColor',[.49 1 .63])

axis([    -0.3516    0.1296   -0.1137    0.1518])
