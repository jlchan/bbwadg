function drawTetWedgeMesh(VX,VY,VZ,EToV)

% WedgeGlobals
K = size(EToV,1);

% plot mesh
hold on
for e = 1:K
    v = EToV(e,:); v = v(v>0); % remove -1 entries -> cruft
    
    if (length(v)==4)
        ids = [1 2 3 1 4 2 3 4];
    else
        ids = [1 2 3 1 4 5 6 4 5 2 3 6];
    end

    text(VX(v)+.1,VY(v)+.1,VZ(v),num2str((1:length(v))'))
    v = v(ids);
    plot3(VX(v),VY(v),VZ(v),'ko-','linewidth',2);    
    
end