function write_LIFTrst

Globals3D
N = 3;
[VX VY VZ] = Nodes3D(N); [VX VY VZ] = xyztorst(VX,VY,VZ);
K = 1; 
EToV = 1:length(VX);

StartUp3D

keyboard
% LIFTr = LIFT*Dr(:,F
% writeMat(LIFTr,'r')
% writeMat(LIFTs,'s')
% writeMat(LIFTt,'t')

function writeMat(LIFT,strname)
fprintf('double p_LIFT%s[%d][%d] = {', strname, Np, Nfp*Nfaces);
for m=1:Np
    fprintf(fid, '{%17.15g ', LIFT(m,1));
    for n=2:Nfp*Nfaces
        fprintf(fid, ', %17.15g ', LIFT(m,n));
    end
    if(m<Np)
        fprintf(fid, '},\n')
    else
        fprintf(fid, '}};\n')
    end
end

