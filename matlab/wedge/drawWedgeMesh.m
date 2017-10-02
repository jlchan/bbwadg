function drawWedgeMesh(VX,VY,VZ,EToV)

% WedgeGlobals
K = size(EToV,1);

% plot mesh
ids = [1 2 3 1 4 5 6 4 5 2 3 6];
hold on
for e = 1:K
    v = EToV(e,:);
text(VX(v)+.1,VY(v)+.1,VZ(v),num2str((1:length(v))'))
    v = v(ids);
    plot3(VX(v),VY(v),VZ(v),'ko-','linewidth',2);
    hold on
    %     for f = 1:Nfaces
    %         enbr = EToE(e,f);
    %         if enbr==e
    %             % do nothing
    %         else
    %             v = EToV(enbr,:);
    %             text(VX(v)+.1,VY(v)+.1,VZ(v),num2str((1:length(v))'))
    %             v = v(ids);
    %             plot3(VX(v),VY(v),VZ(v),'ko-','linewidth',2);
    %         end
    %     end
    %     pause
end