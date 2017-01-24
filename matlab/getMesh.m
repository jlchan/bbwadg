function mesh = getMesh(mesh)

Globals2D
global xq yq

if nargin==0
    mesh.VX = VX;
    mesh.VY = VY;
    mesh.x = x;
    mesh.y = y;
    mesh.EToV = EToV;
    mesh.EToE = EToE;
    mesh.EToF = EToF;
    mesh.K = K;
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.Fscale = Fscale;
    mesh.J = J;
    mesh.sJ = sJ;    
    mesh.mapM = mapM;
    mesh.mapP = mapP;
    mesh.mapB = mapB;
    mesh.vmapM = vmapM;
    mesh.vmapP = vmapP;
    mesh.vmapB = vmapB;
    mesh.xq = xq;
    mesh.yq = yq;
else
    VX = mesh.VX;
    VY = mesh.VY;
    x = mesh.x;
    y = mesh.y;
    EToV = mesh.EToV;
    EToE = mesh.EToE;
    EToF = mesh.EToF;
    K = mesh.K;
    nx = mesh.nx;
    ny = mesh.ny;
    Fscale = mesh.Fscale;
    J = mesh.J;
    sJ = mesh.sJ;    
    mapM = mesh.mapM;
    mapP = mesh.mapP;
    mapB = mesh.mapB;
    vmapM = mesh.vmapM;
    vmapP = mesh.vmapP;
    vmapB = mesh.vmapB;

    xq = mesh.xq;
    yq = mesh.yq;

end

