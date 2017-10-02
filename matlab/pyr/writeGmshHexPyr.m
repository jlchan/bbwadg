% write hex-pyr mesh to GMSH

for K1D = [26];
    [VX VY VZ EToV EToE EToF] = makeHexPyrMesh(K1D);
    
    ids = abs(abs(VX)-1)>1e-2 & abs(abs(VY)-1)>1e-2 & abs(abs(VZ)-1)>1e-2;
    a = .0/K1D;
    VX(ids) = VX(ids) + a*randn(size(VX(ids)));
    VY(ids) = VY(ids) + a*randn(size(VX(ids)));
    VZ(ids) = VZ(ids) + a*randn(size(VX(ids)));
    
    VX = (VX+1)/2;    VY = (VY+1)/2;    VZ = (VZ+1)/2;
    K = size(EToV,1);
    
    
    %fID =
    %fopen(sprintf('../cppV2/benchmarks/elem_timings/pyr/hex_pyr%d.msh',K1D),'w');
    fID = fopen(sprintf('../../cppV2/benchmarks/Cube/hex_pyr_aff%d.msh',K1D),'w');
    
    fprintf(fID,'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
    fprintf(fID,'$Nodes\n')
    fprintf(fID,sprintf('%d\n',length(VX)));
    for i = 1:length(VX)
        fprintf(fID,sprintf('%d %f %f %f\n',i,VX(i),VY(i),VZ(i)));
    end
    fprintf(fID,'$EndNodes\n')
    fprintf(fID,'$Elements\n%d\n',K)
    for e = 1:K
        %fprintf(fID,sprintf('%d 7 2 200 20 %d %d %d %d %d\n',e,EToV(e,1),EToV(e,2),EToV(e,3),EToV(e,4),EToV(e,5)));
        fprintf(fID,sprintf('%d 7 2 0 1 %d %d %d %d %d\n',e,EToV(e,1),EToV(e,2),EToV(e,3),EToV(e,4),EToV(e,5)));
    end
    fprintf(fID,'$EndElements\n')
    fclose(fID);
end