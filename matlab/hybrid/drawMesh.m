function drawMesh(Kvec)

hybridgGlobals3D

if nargin==0
    Kvec = 1:K;
end
hold on
for e = Kvec %1:K
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    switch NvK
        case 4 % tet
            fids = fvT;
        case 5 % pyr
            fids = fvP;            
        case 6 % wedge            
            fids = fvW;
        case 8 % hex            
            fids = fvH;
    end
    Nfaces = length(fids);
    for f = 1:Nfaces
        ids = [fids{f} fids{f}(1)];                
        plot3(VX(EToV(e,ids)),VY(EToV(e,ids)),VZ(EToV(e,ids)),'k-','linewidth',2)
    end
end