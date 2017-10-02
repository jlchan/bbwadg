% write hex-pyr mesh to GMSH

for K1D = 3;
    [VX VY VZ EToV] = makeHexWedgeMesh(K1D);
    K = size(EToV,1);
    
    fID = fopen(sprintf('../../cppV2/benchmarks/elem_timings/prism/hex_pri%d.msh',K1D),'w');
    
    fprintf(fID,'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
    fprintf(fID,'$Nodes\n')
    fprintf(fID,sprintf('%d\n',length(VX)));
    for i = 1:length(VX)
        fprintf(fID,sprintf('%d %f %f %f\n',i,VX(i),VY(i),VZ(i)));
    end
    fprintf(fID,'$EndNodes\n')
    fprintf(fID,'$Elements\n%d\n',K)
    for e = 1:K
        % 6 for pri
        fprintf(fID,sprintf('%d 6 2 200 20 %d %d %d %d %d %d\n',e,EToV(e,1),EToV(e,2),EToV(e,3),EToV(e,4),EToV(e,5),EToV(e,6)));
    end
    fprintf(fID,'$EndElements\n')
    fclose(fID);
end