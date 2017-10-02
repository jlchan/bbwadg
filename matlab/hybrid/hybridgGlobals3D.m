% Purpose: declare pyramid global variables

global N K Nv
global VX VY VZ
global EToE EToF EToV

% VDMs and derivs for volume 
global x y z % cubature nodes
global xf yf zf wf  % cell of Vfq{1:5} - faces

global NpMax NcMax NfcMax

global NpH NcH NfcH
global VH VHr VHs VHt VHf 
global hexK KH xqH yqH zqH 
global fvH fidsH

global NpW NcW NfcW
global VW VWr VWs VWt VWf 
global wedgK KW xqW yqW zqW 
global wW % for wedge LSC-DG
global fvW fidsW
global MW % wedge full mass matrices

global NpP NcP NfcP
global VP VPr VPs VPt VPf 
global pyrK KP xqP yqP zqP 
global fvP fidsP

global NpT NcT NfcT
global VT VTr VTs VTt VTf 
global tetK KT xqT yqT zqT 
global fvT fidsT

% for nodal tets
global xT yT zT rxT sxT txT ryT syT tyT rzT szT tzT JTet
global VTnodal VTq MThat invMThat % quadrature VDMs/mass matrix for nodal tet
% for nodal wedges and tets
global SmaskW SmaskT

% % for plotting
% global xp yp zp
% global VpH VpW VpP VpT sqJp
% global NplotH NplotW NplotP NplotT 

global nx ny nz Fscale sJ  
global mapM mapP mapB

global wJ wsJ invM LIFTT LIFTW xfqW yfqW zfqW 

% for LSC-DG
global J Jsurf Jr Js Jt 
global JW JWf sqJW sqJWf JWr JWs JWt % wedge

% for estimating CFL
global JsB % sJ = has reference face scalings, JsB = no ref face scaling

% for GMSH higher order visualization. 
global rpH spH tpH FH PH VDMH
global rpW spW tpW FW PW VDMW
global rpP spP tpP FP PP VDMP
global rpT spT tpT FT PT VDMT








global rx ry rz sx sy sz tx ty tz
global rk4a rk4b rk4c NODETOL

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];
		 