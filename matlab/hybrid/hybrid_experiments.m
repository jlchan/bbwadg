clear

hybridgGlobals3D

hybrid_mesh
% hex_mesh
% prism_mesh

K = size(EToV,1);
N = 1;

hybridgStartUp

% % make periodic bcs - match node positions
% inIds = abs(xf+1)<NODETOL; outIds = abs(xf-1)<NODETOL;
% 
% [yMp yPp] = meshgrid(yf(inIds),yf(outIds));
% [zMp zPp] = meshgrid(zf(inIds),zf(outIds));
% D = (yMp-yPp').^2 + (zMp-zPp').^2; 
% [id1, id2] = find(abs(D)<NODETOL);
% 
% inIds = inIds(id1); outIds = outIds(id2);
% mapP(inIds) = mapM(outIds); mapP(outIds) = mapM(inIds);

f = @(x,y,z) x + 2*y + 3*z;

%% test L2 projection of polynomials

bT = VT'*(wJT.*f(xT,yT,zT));
bW = VW'*(wJW.*f(xW,yW,zW));
bP = VP'*(wJP.*f(xP,yP,zP));
bH = VH'*(wJH.*f(xH,yH,zH));

u = zeros(NpMax,K);
uT  = invM(1:NpT,tetK).*bT;
uP  = invM(1:NpP,pyrK).*bP;
uW  = invM(1:NpW,wedgK).*bW;
uH  = invM(1:NpH,hexK).*bH;

disp(sprintf('hex error = %f',  norm(VH*uH - f(xH,yH,zH))))
disp(sprintf('wedge error = %f',norm(VW*uW - f(xW,yW,zW))))
disp(sprintf('pyr error = %f',  norm(VP*uP - f(xP,yP,zP))))
disp(sprintf('tet error = %f',  norm(VT*uT - f(xT,yT,zT))))

% hex
rxH = rx(1:NcH,hexK); sxH = sx(1:NcH,hexK); txH = tx(1:NcH,hexK);
ryH = ry(1:NcH,hexK); syH = sy(1:NcH,hexK); tyH = ty(1:NcH,hexK);
rzH = rz(1:NcH,hexK); szH = sz(1:NcH,hexK); tzH = tz(1:NcH,hexK);

% wedge
rxW = rx(1:NcW,wedgK); sxW = sx(1:NcW,wedgK); txW = tx(1:NcW,wedgK);
ryW = ry(1:NcW,wedgK); syW = sy(1:NcW,wedgK); tyW = ty(1:NcW,wedgK);
rzW = rz(1:NcW,wedgK); szW = sz(1:NcW,wedgK); tzW = tz(1:NcW,wedgK);

% pyr
rxP = rx(1:NcP,pyrK); sxP = sx(1:NcP,pyrK); txP = tx(1:NcP,pyrK);
ryP = ry(1:NcP,pyrK); syP = sy(1:NcP,pyrK); tyP = ty(1:NcP,pyrK);
rzP = rz(1:NcP,pyrK); szP = sz(1:NcP,pyrK); tzP = tz(1:NcP,pyrK);

% tet
rxT = rx(1:NcT,tetK); sxT = sx(1:NcT,tetK); txT = tx(1:NcT,tetK);
ryT = ry(1:NcT,tetK); syT = sy(1:NcT,tetK); tyT = ty(1:NcT,tetK);
rzT = rz(1:NcT,tetK); szT = sz(1:NcT,tetK); tzT = tz(1:NcT,tetK);

uHx = rxH.*(VHr*uH) + sxH.*(VHs*uH) + txH.*(VHt*uH);
uHy = ryH.*(VHr*uH) + syH.*(VHs*uH) + tyH.*(VHt*uH);
uHz = rzH.*(VHr*uH) + szH.*(VHs*uH) + tzH.*(VHt*uH);
disp(sprintf('dudx hex error = %f',norm(uHx - 1)))
disp(sprintf('dudy hex error = %f',norm(uHy - 2)))
disp(sprintf('dudz hex error = %f',norm(uHz - 3)))

uWx = rxW.*(VWr*uW) + sxW.*(VWs*uW) + txW.*(VWt*uW);
uWy = ryW.*(VWr*uW) + syW.*(VWs*uW) + tyW.*(VWt*uW);
uWz = rzW.*(VWr*uW) + szW.*(VWs*uW) + tzW.*(VWt*uW);
disp(sprintf('dudx wedge error = %f',norm(uWx - 1)))
disp(sprintf('dudy wedge error = %f',norm(uWy - 2)))
disp(sprintf('dudz wedge error = %f',norm(uWz - 3)))

uPx = rxP.*(VPr*uP) + sxP.*(VPs*uP) + txP.*(VPt*uP);
uPy = ryP.*(VPr*uP) + syP.*(VPs*uP) + tyP.*(VPt*uP);
uPz = rzP.*(VPr*uP) + szP.*(VPs*uP) + tzP.*(VPt*uP);
disp(sprintf('dudx pyr error = %f',norm(uPx - 1)))
disp(sprintf('dudy pyr error = %f',norm(uPy - 2)))
disp(sprintf('dudz pyr error = %f',norm(uPz - 3)))

uTx = rxT.*(VTr*uT) + sxT.*(VTs*uT) + txT.*(VTt*uT);
uTy = ryT.*(VTr*uT) + syT.*(VTs*uT) + tyT.*(VTt*uT);
uTz = rzT.*(VTr*uT) + szT.*(VTs*uT) + tzT.*(VTt*uT);
disp(sprintf('dudx tet error = %f',norm(uTx - 1)))
disp(sprintf('dudy tet error = %f',norm(uTy - 2)))
disp(sprintf('dudz tet error = %f',norm(uTz - 3)))

%% check jumps/fluxes

% get element boundary nodes
u_surface = zeros(NfcMax,K);
u_surface(1:NfcH,hexK ) = VHf*uH;
u_surface(1:NfcW,wedgK) = VWf*uW;
u_surface(1:NfcP,pyrK ) = VPf*uP;
u_surface(1:NfcT,tetK ) = VTf*uT;

u_diff = zeros(NfcMax,K);
u_diff(1:NfcH,hexK)  = u_surface(mapM(1:NfcH,hexK)) -u_surface(mapP(1:NfcH,hexK));
u_diff(1:NfcW,wedgK) = u_surface(mapM(1:NfcW,wedgK))-u_surface(mapP(1:NfcW,wedgK));
u_diff(1:NfcP,pyrK)  = u_surface(mapM(1:NfcP,pyrK)) -u_surface(mapP(1:NfcP,pyrK));
u_diff(1:NfcT,tetK)  = u_surface(mapM(1:NfcT,tetK)) -u_surface(mapP(1:NfcT,tetK));
norm(u_diff,'fro')
% uTf = VTf*uT;
% uTf = VTf*uT;
% uTf = VTf*uT;
color_line3(xf,yf,zf,u_diff,'.');
colorbar

% nx(mapB) 