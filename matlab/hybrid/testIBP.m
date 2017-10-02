clear

hybridgGlobals3D

% hybrid_mesh
% hex_mesh
% prism_mesh

r1 = [-1  1  1 -1 -1  1  1 -1]'; s1 = [-1 -1  1  1 -1 -1  1  1]';
t1 = [-1 -1 -1 -1  1  1  1  1]';
rstH = [r1 s1 t1];

% wedge face vertex ordering
u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
r1 = v(:); s1 = w(:); t1 = u(:); % flipping coordinates for Gmsh
rstW = [r1 s1 t1];

% pyramid face vertex ordering
r1 = [ -1   1   1  -1  -1 ]';s1 = [ -1  -1   1   1  -1 ]';
t1 = [ -1  -1  -1  -1   1 ]';
rstP = [r1 s1 t1];

% tet face vertex orderings
r1 = [-1  1 -1 -1 ]'; s1 = [-1 -1  1 -1 ]'; t1 = [-1 -1 -1  1 ]';
rstT = [r1 s1 t1];

% choose type of element
rst = rstT; 

VX = rst(:,1); VY = rst(:,2); VZ = rst(:,3);
EToV = (1:length(VX));
EToE = 0*EToV;
EToF = 0*EToV;
K = size(EToV,1);
N = 2;

hybridgStartUp

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

wJH = wJ(1:NcH,hexK);  wsJH = wsJ(1:NfcH,hexK);  nxH = nx(1:NfcH,hexK);
wJW = wJ(1:NcW,wedgK); wsJW = wsJ(1:NfcW,wedgK); nxW = nx(1:NfcW,wedgK);
wJP = wJ(1:NcP,pyrK);  wsJP = wsJ(1:NfcP,pyrK);  nxP = nx(1:NfcP,pyrK);
wJT = wJ(1:NcT,tetK);  wsJT = wsJ(1:NfcT,tetK);  nxT = nx(1:NfcT,tetK);

% =============== projection test ============

f = @(x,y,z) x + y + z;

uH  = invM(1:NpH,hexK).*(VH'*(wJH.*f(xH,yH,zH)));
uW  = invM(1:NpW,wedgK).*(VW'*(wJW.*f(xW,yW,zW)));
uP  = invM(1:NpP,pyrK).*(VP'*(wJP.*f(xP,yP,zP)));
uT  = invM(1:NpT,tetK).*(VT'*(wJT.*f(xT,yT,zT)));

if Nv==4
    
    rx = rxT; sx = sxT; tx = txT; nx = nxT;
    wJ = wJT; wsJ = wsJT; 
    V = VT; Vr = VTr; Vs = VTs; Vt = VTt; Vf = VTf;
    u = uT;
    
elseif Nv==5
    
    rx = rxP; sx = sxP; tx = txP; nx = nxP;
    wJ = wJP; wsJ = wsJP;    
    V = VP; Vr = VPr; Vs = VPs; Vt = VPt; Vf = VPf;    
    u = uP;
    
elseif Nv==6
    
    rx = rxW; sx = sxW; tx = txW; nx = nxW;
    wJ = wJW; wsJ = wsJW;
    V = VW; Vr = VWr; Vs = VWs; Vt = VWt; Vf = VWf;
    u = uW;
    
elseif Nv==8
    
    rx = rxH; sx = sxH; tx = txH; nx = nxH;
    wJ = wJH; wsJ = wsJH;
    V = VH; Vr = VHr; Vs = VHs; Vt = VHt;  Vf = VHf;
    u = uH;
    
end

% integral of v * dudx
dudx = rx.*(Vr*u) + sx.*(Vs*u) + tx.*(Vt*u);
S1 = -V'*(wJ.*dudx);
% integral of (-dvdx * u)  +  surf_int (v * u)
% dvdx = dvdr * rx * u + dvds * sx * u + dvdt * tx * u

uc = wJ.*(V*u);  
S2 = (Vr'*(rx.*uc) + Vs'*(sx.*uc) + Vt'*(tx.*uc)) + ...
    - Vf'*(nx.*wsJ.*(Vf*u));

norm(S1-S2,'fro')
