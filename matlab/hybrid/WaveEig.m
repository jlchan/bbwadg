function WaveEig

hybridgGlobals3D
hybridgGlobalFlags
useSEM = 1; useLSC = 1; useNodalTets = 1; useSkewRHS = 0;

% prism_mesh
% hybrid_mesh
% EToV(EToV==0) = -1; EToV = EToV'; EToE = EToE'; EToF = EToF'; % if
% hybrid_mesh2
% pyr_mesh
% prism_mesh2

% hex2
VX = [    -1     1    1     -1    -1     1    1     -1]'; VY = [    -1    -1     1     1    -1    -1     1     1]';  VZ = [    -1    -1    -1    -1     1     1     1     1]';

ids = abs(abs(VX)-1)>1e-2 & abs(abs(VY)-1)>1e-2 & abs(abs(VZ)-1)>1e-2;
a = .2;
VX(ids) = VX(ids) + a*randn(size(VX(ids))); 
VY(ids) = VY(ids) + a*randn(size(VX(ids))); 
VZ(ids) = VZ(ids) + a*randn(size(VX(ids)));

% EToV = 1:8; EToE = 0*EToV; EToF = 0*EToV;

K = size(EToV,1);

% drawMesh
% % text(VX,VY,VZ,num2str((1:length(VX))'))
% axis equal
% print(gcf,'-depsc','../docs/figs/hybridMesh.eps')
% keyboard
Nmax = 4;
Nvec = [Nmax];
Cmat = hsv(Nmax);
spectra = cell(Nmax,1);
spectra_sym = cell(Nmax,1);
spectra_skew = cell(Nmax,1);

ARHS = cell(Nmax,1);
for N = Nvec
    
    hybridgStartUp
    
    [A M Minv] = buildOperator();
    ARHS{N} = A; % M\A 
            
    MA = M*A; % returns back just A
    Asym = .5*Minv*(MA + MA');
    Askew = .5*Minv*(MA - MA');
    
    lam = eig(A);
    lam_sym = eig(Asym);
    lam_skew = eig(Askew);
    
    spectra{N} = lam;
    spectra_sym{N} = lam_sym;
    spectra_skew{N} = lam_skew;
       
    N
end

keyboard

% %%
% rxK = rx(1:NcH,hexK); ryK = ry(1:NcH,hexK); rzK = rz(1:NcH,hexK);
% sxK = sx(1:NcH,hexK); syK = sy(1:NcH,hexK); szK = sz(1:NcH,hexK);
% txK = tx(1:NcH,hexK); tyK = ty(1:NcH,hexK); tzK = tz(1:NcH,hexK);
% JK = J(1:NcH,hexK); JsK = JsB(1:NfcH,hexK);
% 
% rm = max([max(abs(rxK)),max(abs(ryK)),max(abs(rzK))]);
% sm = max([max(abs(sxK)),max(abs(syK)),max(abs(szK))]);
% tm = max([max(abs(txK)),max(abs(tyK)),max(abs(tzK))]);
% max([rm,sm,tm]) * max(JK)*min(JK)
% max(JsK) 
%% compute trace constant bounds
rho = zeros(Nmax,1); rhoSym = zeros(Nmax,1); rhoSkew = zeros(Nmax,1);
bC = zeros(Nmax,1); C = zeros(Nmax,1);
biC = zeros(Nmax,1);

% precomputed trace constants
PriConst = [9.926135933275306  18.563576705381948  29.030325215439685  42.988345972840094 58.802145509223905  78.006158337859318  99.271490513773074];
PyrConst = 1.0e+02 *[ 0.116837586022311   0.208869058826159   0.328373010378272   0.475870076602500    0.651676634297143   0.855997828354466   1.088964369062309];
TetConst = [ 12.222306675423649  20.461046774606864  29.175757913298817  41.653800010846517 54.447306071124011  71.104042033768295  88.317462093087173];

% precomputed markov constants
CMH =[    9.0000   22.4375   27.2336   30.2598   32.1556   33.5375   34.6033];
CMW = 1.0e+03 * [   0.012000000000000   0.054270509831248   0.142634745966637   0.308340065156952   0.585888670274778   1.021635581744581   1.663852290534619];
CMP =   1.0e+03 *[  0.012916666666667   0.060050988028815   0.175509522381149   0.405433613335050 0.809951819422284   1.460934071720827   2.442261122895831];
CMT = 1.0e+03 * [0.020000000000000   0.078616173517937   0.195578485633278   0.403911155894613 0.744853548366754   1.265543961520880   2.021087045974716];

for N = 1:Nmax
    
    hybridgStartUp

    %% hex    
    CNH = 3*(N+1).*(N+2)/2; % trace constant over the surface  
    bCH = 0; biCH = 0;
    if ~isempty(hexK)
        JsH = JsB(1:NfcH,hexK); JH = J(1:NcH,hexK);
        bCH = CNH * max(JsH(:)) * max(1./JH(:));        
            
        rxK = rx(hexK); ryK = ry(hexK); rzK = rz(hexK);
        sxK = sx(hexK); syK = sy(hexK); szK = sz(hexK);
        txK = tx(hexK); tyK = ty(hexK); tzK = tz(hexK);
        Crst = max([max(rxK(:)),max(ryK(:)),max(rzK(:)),...
            max(sxK(:)),max(syK(:)),max(szK(:)),...
            max(txK(:)),max(tyK(:)),max(tzK(:))]);
        biCH = (Crst*sqrt(CMH(N)*max(JH(:))/min(JH(:))) + CNH * max(JsH(:))) * max(1./JH(:));

    end    

    CH = 0;
    for ee = 1:length(hexK)
        e = hexK(ee);
        Mf = VHf'*spdiag(wsJ(1:NfcH,e))*VHf; % incorporate face scalings
        M = VH'*spdiag(wJ(1:NcH,e))*VH;
        CH = max(CH,max(abs(eig(M\Mf))));
    end
    
    %% wedge    
    CNW = PriConst(N); 
%     CNW = (sqrt(2)+3)*(N+1)*(N+2); %(N+1).*(N+2)/2;
    bCW = 0; biCW = 0;
    if ~isempty(wedgK)
        JsW = JsB(1:NfcW,wedgK);  JW = J(1:NcW,wedgK);
        bCW = CNW * max(JsW(:))*max(1./JW(:));
        
        rxK = rx(wedgK); ryK = ry(wedgK); rzK = rz(wedgK);
        sxK = sx(wedgK); syK = sy(wedgK); szK = sz(wedgK);
        txK = tx(wedgK); tyK = ty(wedgK); tzK = tz(wedgK);
        Crst = max([max(rxK(:)),max(ryK(:)),max(rzK(:)),...
            max(sxK(:)),max(syK(:)),max(szK(:)),...
            max(txK(:)),max(tyK(:)),max(tzK(:))]);
        biCW = (Crst*sqrt(CMW(N)*max(JW(:))/min(JW(:))) + CNW * max(JsW(:))) * max(1./JW(:));
    end    
    
    CW = 0;
    for ee = 1:length(wedgK)
        e = wedgK(ee);
        Mf = VWf'*spdiag(wsJ(1:NfcW,e))*VWf; % incorporate face scalings
        M = VW'*spdiag(wJ(1:NcW,e))*VW;
        CW = max(CW,max(abs(eig(M\Mf))));
    end
    
    %% pyr
    CNP = PyrConst(N); 
%     CNP = (sqrt(2)+2)*(N+1)*(N+3); 
    bCP = 0; biCP = 0;
    if ~isempty(pyrK)
        JsP = JsB(1:NfcP,pyrK); JP = J(1:NcP,pyrK);    
        bCP = CNP * max(JsP(:)) * max(1./JP(:));
        
        rxK = rx(pyrK); ryK = ry(pyrK); rzK = rz(pyrK);
        sxK = sx(pyrK); syK = sy(pyrK); szK = sz(pyrK);
        txK = tx(pyrK); tyK = ty(pyrK); tzK = tz(pyrK);
        Crst = max([max(rxK(:)),max(ryK(:)),max(rzK(:)),...
            max(sxK(:)),max(syK(:)),max(szK(:)),...
            max(txK(:)),max(tyK(:)),max(tzK(:))]);
        biCP = (Crst*sqrt(CMP(N)*max(JP(:))/min(JP(:))) + CNP * max(JsP(:))) * max(1./JP(:));
    end
    
    CP = 0;
    for ee = 1:length(pyrK)
        e = pyrK(ee);
        Mf = VPf'*spdiag(wsJ(1:NfcP,e))*VPf; % incorporate face scalings
        M = VP'*spdiag(wJ(1:NcP,e))*VP;
        CP = max(CP,max(abs(eig(M\Mf))));
    end
    
    %% tet
    
    CNT = TetConst(N); 
%     CNT = (sqrt(3)+3)*(N+1).*(N+3)/2;
    bCT = 0; biCT = 0;
    if ~isempty(tetK)
    JsT = JsB(1:NfcT,tetK); JT = J(1:NcT,tetK);    
        bCT = CNT * max(JsT(:)) * max(1./JT(:));
        
        rxK = rx(tetK); ryK = ry(tetK); rzK = rz(tetK);
        sxK = sx(tetK); syK = sy(tetK); szK = sz(tetK);
        txK = tx(tetK); tyK = ty(tetK); tzK = tz(tetK);
        Crst = max([max(rxK(:)),max(ryK(:)),max(rzK(:)),...
            max(sxK(:)),max(syK(:)),max(szK(:)),...
            max(txK(:)),max(tyK(:)),max(tzK(:))]);
        biCT = (Crst*sqrt(CMT(N)*max(JT(:))/min(JT(:))) + CNP * max(JsT(:))) * max(1./JT(:));
    end
    
    CT = 0;
    for ee = 1:length(tetK)
        e = tetK(ee);
        Mf = VTf'*spdiag(wsJ(1:NfcT,e))*VTf; % incorporate face scalings
        M = VT'*spdiag(wJ(1:NcT,e))*VT;
        CT = max(CT,max(abs(eig(M\Mf))));
    end
    
    
    %% true/estimated trace constants
    rho(N) = max(abs(spectra{N}));
    rhoSym(N) = max(abs(spectra_sym{N}));    
    rhoSkew(N) = max(abs(spectra_skew{N}));    
    bC(N) = max([bCH bCW bCP bCT]);
    biC(N)  = max([biCH biCW biCP biCT]);    
    C(N)  = max([CH CW CP CT]);    
    
    fprintf('done with N = %i\n',N);
end


%%
figure
hold on
for N = Nmax:-1:1
    lam = spectra{N};
    lam_sym = spectra_sym{N};
    lam_skew = spectra_skew{N};
        
    % bound based on skew/sym parts    
    rlam = linspace(min(real(lam_sym)),max(real(lam_sym)),50);
    ilam = linspace(min(imag(lam_skew)),max(imag(lam_skew)),50);
    box1 = [rlam + 1i*max(ilam),nan, rlam + 1i*min(ilam)];
    box2 = [1i*ilam + max(rlam),nan, 1i*ilam + min(rlam)];
    boxLegend(N) = plot([box1 box2],'--','color',Cmat(N,:),'linewidth',2,'DisplayName','Sym/skew bounds');
        
    eigLegend(N) = plot(lam,'.','color',Cmat(N,:),'markersize',16,'DisplayName',sprintf('Spectra for N = %d',N));
    
%     % ==== box bounds
%     % bound based on exact trace ineq    
%     ilam = linspace(-biC(N),biC(N),10);
%     rlam = linspace(-bC(N),0,10);    
%     box1 = [rlam + 1i*max(ilam),nan, rlam + 1i*min(ilam)];
%     box2 = [1i*ilam + max(rlam),nan, 1i*ilam + min(rlam)];
%     plot([box1 box2],'s--','color',Cmat(N,:),'linewidth',2,'DisplayName','Derived bounds')
    
end

legend(eigLegend)
axis equal
axis([  -55.5873    1.6144  -22.2046   22.9108])
set(gca,'fontsize',14)

% zeroAx = 0 + 1i*ilam;
% plot(zeroAx,'k-.','linewidth',2,'DisplayName','Imaginary Axis')

% print(gcf,'-depsc','../docs/figs/symSkewEigBounds.eps')
%%
keyboard

print(gcf,'-depsc','../docs/figs/symSkewEigHybrid.eps')
return


function [A M Minv] = buildOperator

hybridgGlobals3D
hybridgGlobalFlags

% get eigenvalues
e1 = zeros(NpMax,K); e2 = zeros(NpMax,K);
e3 = zeros(NpMax,K); e4 = zeros(NpMax,K);

% pick out nonzero entries of dof matrix
e1(1:NpH,hexK) = 1; e1(1:NpW,wedgK) = 1;
e1(1:NpP,pyrK) = 1; e1(1:NpT,tetK) = 1;
inds = find(e1);

Ndofs = NpH*KH + NpW*KW + NpP*KP + NpT*KT;
A = zeros(4*Ndofs);
for i = 1:4*Ndofs
    e = zeros(4*Ndofs,1);     e(i) = 1;
    e1v = e(1:Ndofs);
    e2v = e(Ndofs + (1:Ndofs));
    e3v = e(2*Ndofs + (1:Ndofs));
    e4v = e(3*Ndofs + (1:Ndofs));
    
    e1(inds) = e1v; e2(inds) = e2v; e3(inds) = e3v; e4(inds) = e4v;
    
    [e1s e2s e3s e4s] = SurfaceInterp(e1,e2,e3,e4);
    
    [rhsp rhsu rhsv rhsw] = RHSall(e1,e2,e3,e4,e1s,e2s,e3s,e4s);
    
    A(0*Ndofs + (1:Ndofs),i) = rhsp(inds);
    A(1*Ndofs + (1:Ndofs),i) = rhsu(inds);
    A(2*Ndofs + (1:Ndofs),i) = rhsv(inds);
    A(3*Ndofs + (1:Ndofs),i) = rhsw(inds);
end
A(abs(A)<1e-8)=0; % A = M\A

% use fact that elements are ordered hex, wedg, pyr, tet
Md = []; iMd = [];
for e = 1:K
    NvK = nnz(EToV(e,:)>0);
    if NvK==4
        Np = NpT;
    elseif NvK==5
        Np=NpP;
    elseif NvK==6
        Np=NpW;
    elseif NvK==8
        Np = NpH;
    end
    Md = [Md; 1./invM(1:Np,e)];
    iMd = [iMd; invM(1:Np,e)];
end
M = diag(Md); M = blkdiag(M,M,M,M); % four times for p,u,v,w
Minv = diag(iMd); Minv = blkdiag(Minv,Minv,Minv,Minv); % four times for p,u,v,w

return





% computes RHS for all element types 
function [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws)

hybridgGlobals3D
hybridgGlobalFlags

% hex
if useSkewRHS
    [rhspK rhsuK rhsvK rhswK] = SkewRHS(p,u,v,w,ps,us,vs,ws,VH,VHr,VHs,VHt,VHf,hexK );
else
    [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VH,VHr,VHs,VHt,VHf,hexK );
end
rhsp(1:NpH,hexK) = rhspK; rhsu(1:NpH,hexK) = rhsuK;
rhsv(1:NpH,hexK) = rhsvK; rhsw(1:NpH,hexK) = rhswK;

% wedges
if useLSC
    [rhspK rhsuK rhsvK rhswK] = LSC_RHS(p,u,v,w,ps,us,vs,ws,...
        VW,VWr,VWs,VWt,VWf,wedgK,sqJW, sqJWf, JWr, JWs, JWt);
else
    [rhspK rhsuK rhsvK rhswK] = Wedge_RHS(p,u,v,w,ps,us,vs,ws,VW,VWr,VWs,VWt,VWf,wedgK);
end
rhsp(1:NpW,wedgK) = rhspK; rhsu(1:NpW,wedgK) = rhsuK;
rhsv(1:NpW,wedgK) = rhsvK; rhsw(1:NpW,wedgK) = rhswK;

% pyrs - always use skew form, same cost but stable if useSEM=1
% if useSkewRHS
    [rhspK rhsuK rhsvK rhswK] = SkewRHS(p,u,v,w,ps,us,vs,ws,VP,VPr,VPs,VPt,VPf,pyrK );
% else
%     [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VP,VPr,VPs,VPt,VPf,pyrK );
% end
rhsp(1:NpP,pyrK ) = rhspK; rhsu(1:NpP,pyrK ) = rhsuK;
rhsv(1:NpP,pyrK ) = rhsvK; rhsw(1:NpP,pyrK ) = rhswK;

% tet - switch w/nodal DG
if useNodalTets
    [rhspK rhsuK rhsvK rhswK] = NodalRHS(p,u,v,w,ps,us,vs,ws,VTr,VTs,VTt,VTf);
else
    [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VT,VTr,VTs,VTt,VTf,tetK );
end
rhsp(1:NpT,tetK ) = rhspK; rhsu(1:NpT,tetK ) = rhsuK;
rhsv(1:NpT,tetK ) = rhsvK; rhsw(1:NpT,tetK ) = rhswK;

% nodal RHS for tets: "collocation" DG
function [rhsp rhsu rhsv rhsw] = NodalRHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    Vr,Vs,Vt,Vf)

hybridgGlobals3D;

if nnz(tetK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Np = NpT; Nfc = size(Vf,1);

% face cubature
wsJK = wsJ(1:Nfc,tetK); 
nxK = nx(1:Nfc,tetK); nyK = ny(1:Nfc,tetK); nzK = nz(1:Nfc,tetK);

% nodal dofs
pK = p(1:Np,tetK); uK = u(1:Np,tetK); vK = v(1:Np,tetK); wK = w(1:Np,tetK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,tetK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxT.*(Vr*uK) + sxT.*(Vs*uK) + txT.*(Vt*uK);
dvdy = ryT.*(Vr*vK) + syT.*(Vs*vK) + tyT.*(Vt*vK);
dwdz = rzT.*(Vr*wK) + szT.*(Vs*wK) + tzT.*(Vt*wK);
divU = dudx + dvdy + dwdz;
rhsp = -divU + invMThat*((Vf'*(wsJK.*flux_p))./JTet);

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxT.*pr + sxT.*ps + txT.*pt;
dpdy = ryT.*pr + syT.*ps + tyT.*pt;
dpdz = rzT.*pr + szT.*ps + tzT.*pt;

flux_uc = wsJK.*flux_u;
rhsu = -dpdx + invMThat*((Vf'*(nxK.*flux_uc))./JTet);
rhsv = -dpdy + invMThat*((Vf'*(nyK.*flux_uc))./JTet);
rhsw = -dpdz + invMThat*((Vf'*(nzK.*flux_uc))./JTet);

return

% strong form RHS with quadrature - OK for hex + pyramids + tets. 
% can do better w/nodal tets
function [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
divU = dudx + dvdy + dwdz;
rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

function [rhsp rhsu rhsv rhsw] = SkewRHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
% flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);
flux_p = .5*p_jump - nu_avg;

% integrated by parts pressure equation
uc = wJK.*(V*uK); vc = wJK.*(V*vK); wc = wJK.*(V*wK);
Vxu = Vr'*(rxK.*uc) + Vs'*(sxK.*uc) + Vt'*(txK.*uc);
Vyv = Vr'*(ryK.*vc) + Vs'*(syK.*vc) + Vt'*(tyK.*vc);
Vzw = Vr'*(rzK.*wc) + Vs'*(szK.*wc) + Vt'*(tzK.*wc);
gradvU = Vxu + Vyv + Vzw;
rhsp = invMK.*(gradvU + Vf'*(wsJK.*flux_p));

% % volume derivative operators
% dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
% dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
% dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
% divU = dudx + dvdy + dwdz;
% rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

function [rhsp rhsu rhsv rhsw] = LSC_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK,...
    sqJK, sqJKf, JKr, JKs, JKt)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); % could also just take typeK(1) - same mass matrix
wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*p_jump - nu_avg;
flux_u = .5*(nu_jump - p_jump);

% % integrated by parts pressure equation
% LSC approach - chain rule derivatives
uc = wJK.*(V*uK).*sqJK.*sqJK;
Vrxu = Vr'*(rxK.*uc) - V'*(JKr.*rxK.*uc);
Vsxu = Vs'*(sxK.*uc) - V'*(JKs.*sxK.*uc);
Vtxu = Vt'*(txK.*uc) - V'*(JKt.*txK.*uc);

vc = wJK.*(V*vK).*sqJK.*sqJK;
Vryv = Vr'*(ryK.*vc) - V'*(JKr.*ryK.*vc);
Vsyv = Vs'*(syK.*vc) - V'*(JKs.*syK.*vc);
Vtyv = Vt'*(tyK.*vc) - V'*(JKt.*tyK.*vc);

wc = wJK.*(V*wK).*sqJK.*sqJK;
Vrzw = Vr'*(rzK.*wc) - V'*(JKr.*rzK.*wc);
Vszw = Vs'*(szK.*wc) - V'*(JKs.*szK.*wc);
Vtzw = Vt'*(tzK.*wc) - V'*(JKt.*tzK.*wc);

Vxu = Vrxu + Vsxu + Vtxu;
Vyv = Vryv + Vsyv + Vtyv;
Vzw = Vrzw + Vszw + Vtzw;
gradvU = Vxu + Vyv + Vzw;
rhsp = invMK.*(gradvU + Vf'*(sqJKf.*wsJK.*flux_p)); % invM = ones for LSC-DG 

% volume derivative operators
% pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;

pr = (Vr*pK - (V*pK).*JKr).*sqJK;
ps = (Vs*pK - (V*pK).*JKs).*sqJK;
pt = (Vt*pK - (V*pK).*JKt).*sqJK;

dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = sqJKf.*wsJK.*flux_u;
sqwJK = sqJK.*wJK; 

rhsu = invMK.*(-V'*(sqwJK.*dpdx) + Vf'*(nxK.*flux_uc)); % invM = ones for LSC-DG
rhsv = invMK.*(-V'*(sqwJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(sqwJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

% strong form RHS with full mass inversion for wedges
function [rhsp rhsu rhsv rhsw] = Wedge_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

% invMK = invM(1:Np,typeK); 
wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
divU = dudx + dvdy + dwdz;
% rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
% rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
% rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
% rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

% full mass matrix inversion
rhsp = (-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));
rhsu = (-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = (-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = (-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));
for e = 1:length(typeK)
    rhse = MW{e}\[rhsp(:,e) rhsu(:,e) rhsv(:,e) rhsw(:,e)];
    rhsp(:,e) = rhse(:,1);    rhsu(:,e) = rhse(:,2);
    rhsv(:,e) = rhse(:,3);    rhsw(:,e) = rhse(:,4);
end

return

% independent of regular DG or LSC-DG
function [p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,...
    Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    p_jump = [];    nu_jump = [];    nu_avg = [];
    return
end

Nfc = size(Vf,1); Np = size(Vf,2);
nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);

mapMK = mapM(1:Nfc,typeK); mapPK = mapP(1:Nfc,typeK);

% flux gather
p_jump  = p_surface(mapPK) - p_surface(mapMK);
u_jump  = u_surface(mapPK) - u_surface(mapMK);
v_jump  = v_surface(mapPK) - v_surface(mapMK);
w_jump  = w_surface(mapPK) - w_surface(mapMK);
nu_jump = nxK.*u_jump + nyK.*v_jump + nzK.*w_jump;

u_avg  = .5*(u_surface(mapPK) + u_surface(mapMK));
v_avg  = .5*(v_surface(mapPK) + v_surface(mapMK));
w_avg  = .5*(w_surface(mapPK) + w_surface(mapMK));
nu_avg = nxK.*u_avg + nyK.*v_avg + nzK.*w_avg;

% apply BCs, todo: try p^+ = - p ^ -, which is energy stable
mapBK = mapMK==mapPK;
p_jump(mapBK) = -2*p_surface(mapMK(mapBK));


function [ps us vs ws] = SurfaceInterp(p,u,v,w)

hybridgGlobals3D
hybridgGlobalFlags

ps = zeros(NfcMax,K);
us = zeros(NfcMax,K);
vs = zeros(NfcMax,K);
ws = zeros(NfcMax,K);

ps(1:NfcH,hexK ) = VHf*p(1:NpH,hexK );
us(1:NfcH,hexK ) = VHf*u(1:NpH,hexK );
vs(1:NfcH,hexK ) = VHf*v(1:NpH,hexK );
ws(1:NfcH,hexK ) = VHf*w(1:NpH,hexK );

if useLSC
    % wedges LSC interp
    ps(1:NfcW,wedgK) = (VWf*p(1:NpW,wedgK)).*sqJWf;
    us(1:NfcW,wedgK) = (VWf*u(1:NpW,wedgK)).*sqJWf;
    vs(1:NfcW,wedgK) = (VWf*v(1:NpW,wedgK)).*sqJWf;
    ws(1:NfcW,wedgK) = (VWf*w(1:NpW,wedgK)).*sqJWf;
else
    ps(1:NfcW,wedgK) = VWf*p(1:NpW,wedgK);
    us(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
    vs(1:NfcW,wedgK) = VWf*v(1:NpW,wedgK);
    ws(1:NfcW,wedgK) = VWf*w(1:NpW,wedgK);    
end

ps(1:NfcP,pyrK ) = VPf*p(1:NpP,pyrK );
us(1:NfcP,pyrK ) = VPf*u(1:NpP,pyrK );
vs(1:NfcP,pyrK ) = VPf*v(1:NpP,pyrK );
ws(1:NfcP,pyrK ) = VPf*w(1:NpP,pyrK );

ps(1:NfcT,tetK ) = VTf*p(1:NpT,tetK );
us(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );
vs(1:NfcT,tetK ) = VTf*v(1:NpT,tetK );
ws(1:NfcT,tetK ) = VTf*w(1:NpT,tetK );

