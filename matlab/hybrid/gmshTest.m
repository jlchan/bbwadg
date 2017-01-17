clear
clearvars -global

hybridgGlobals3D
hybridgGlobalFlags
useSEM = 1; useLSC = 1; useNodalTets = 1; useSkewRHS = 1;

% hybrid_mesh
% prism2; EToV = EToV'; EToE = EToE'; EToF = EToF'; 
% prism4; EToV = EToV'; EToE = EToE'; EToF = EToF'; 
% pyr_mesh6; EToV = EToV'; EToE = EToE'; EToF = EToF'; 

% VX = [ -1   1   1  -1   0]'; 
% VY = [ -1  -1   1   1   0]';  
% VZ = [  0   0   0   0   1]';
% EToV = 1:length(VX); EToE = 0*EToV; EToF = 0*EToV;

prism_Twist2;

% VX = [    -1     1    1     -1    -1     1    1     -1]'; 
% VY = [    -1    -1    1     1    -1    -1     1     1]';  
% VZ = [    -1    -1    -1    -1    1     1     1     1]';
% EToV = 1:length(VX); EToE = 0*EToV; EToF = 0*EToV;

% hex_mesh

% cube_mesh

K = size(EToV,1);
N = 3;

hybridgStartUp
% drawMesh

%% set initial conditions

% exact sol
uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);
f = @(x,y,z) uexf(x,y,z,0);
% f = @(x,y,z) (1-x).*(1+x).*(1-y).*(1+y).*(1-z).*(1+z);

b(1:NpH,hexK ) = VH'*(wJ(1:NcH,hexK ).*f(xqH,yqH,zqH));
if useLSC
    b(1:NpW,wedgK) = VW'*(sqJW.*wJ(1:NcW,wedgK).*f(xqW,yqW,zqW));
else
    b(1:NpW,wedgK) = VW'*(wJ(1:NcW,wedgK).*f(xqW,yqW,zqW));
end
b(1:NpP,pyrK ) = VP'*(wJ(1:NcP,pyrK ).*f(xqP,yqP,zqP));
b(1:NpT,tetK ) = VT'*(wJ(1:NcT,tetK ).*f(xqT,yqT,zqT));

p = invM.*b; % LSC inverse built into invM - extra storage but OK


%% eval at points

xp = zeros(NpMax,K); yp = zeros(NpMax,K); zp = zeros(NpMax,K);
pplot = zeros(NpMax,K);
for e = 1:K
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    switch NvK
        
        case 4 % tet
            
            [xpK,ypK,zpK] = tet_geom_factors(VX(v),VY(v),VZ(v),rpT,spT,tpT);
            pplot(1:length(rpT),e) = VDMT * p(1:NpT,e);
            NpK = length(rpT);
        
        case 5 % pyr
            
            [xpK,ypK,zpK] = pyr_geom_factors(VX(v),VY(v),VZ(v),rpP,spP,tpP);
            pplot(1:length(rpP),e) = VDMP * p(1:NpP,e);
            NpK = length(rpP);
            
        case 6 % wedge
            
            [xpK,ypK,zpK, ~,~,~,~,~,~,~,~,~,JpW] = ...
                wedge_geom_factors(VX(v),VY(v),VZ(v),rpW,spW,tpW);
            
            % LSC-DG wedge.
            pplot(1:length(rpW),e) = (VDMW * p(1:NpW,e))./sqrt(JpW);
            NpK = length(rpW);
            
        case 8 % hex
            
            [xpK,ypK,zpK] = hex_geom_factors(VX(v),VY(v),VZ(v),rpH,spH,tpH);
            pplot(1:length(rpH),e) = VDMH * p(1:NpH,e);
            NpK = length(rpH);
            
    end
    xp(1:NpK,e) = xpK;    yp(1:NpK,e) = ypK;    zp(1:NpK,e) = zpK;
end
color_line3(xp,yp,zp,pplot,'.')
% break

%%
fID = fopen(sprintf('p_output.msh'),'w');

fprintf(fID,'$MeshFormat\n');
fprintf(fID,'2.2 0 8\n');
fprintf(fID,'$EndMeshFormat\n');

fprintf(fID,'$InterpolationScheme\n');
fprintf(fID,'"Nodal_Poly_Interpolation"\n'); %"name"
fprintf(fID,'4\n'); % number-of-element-topologies

% tet = elem 4
fprintf(fID,'5 2\n'); % elm-topology, number-of-interpolation-matrices
fprintf(fID,'%d %d\n',size(FT,1),size(FT,2));%num-rows num-columns value
fprintf(fID, [repmat('%f ', 1, size(FT, 2)) '\n'], FT);
fprintf(fID,'%d %d\n',size(PT,1),size(PT,2));%num-rows num-columns value
fprintf(fID, [repmat('%d ', 1, size(PT, 2)) '\n'], PT');

% pyr = elem 6
fprintf(fID,'6 2\n'); % elm-topology, number-of-interpolation-matrices
fprintf(fID,'%d %d\n',size(FP,1),size(FP,2));%num-rows num-columns value
fprintf(fID, [repmat('%f ', 1, size(FP, 2)) '\n'], FP);
fprintf(fID,'%d %d\n',size(PP,1),size(PP,2));%num-rows num-columns value
fprintf(fID, [repmat('%d ', 1, size(PP, 2)) '\n'], PP');
 
% wedge = elem 7
fprintf(fID,'7 2\n'); % elm-topology, number-of-interpolation-matrices
fprintf(fID,'%d %d\n',size(FW,1),size(FW,2));%num-rows num-columns value
fprintf(fID, [repmat('%f ', 1, size(FW, 2)) '\n'], FW);
fprintf(fID,'%d %d\n',size(PW,1),size(PW,2)); % num-rows num-columns value
fprintf(fID, [repmat('%d ', 1, size(PW, 2)) '\n'], PW');
 
% hex = elem 5
fprintf(fID,'8 2\n'); % elm-topology, number-of-interpolation-matrices
fprintf(fID,'%d %d\n',size(FH,1),size(FH,2)); % num-rows num-columns value
fprintf(fID, [repmat('%f ', 1, size(FH, 2)) '\n'], FH);
fprintf(fID,'%d %d\n',size(PH,1),size(PH,2)); % num-rows num-columns value
fprintf(fID, [repmat('%d ', 1, size(PH, 2)) '\n'], PH');

fprintf(fID,'$EndInterpolationScheme\n');



% hex
fprintf(fID,'$ElementNodeData\n');
fprintf(fID,'2\n"hexData"\n"Nodal_Poly_Interpolation"\n'); % #-string-tags, "string-tag", interp scheme
fprintf(fID,'1\n 0.1\n'); % number-of-real-tags, < real-tag > - fake time tag
fprintf(fID,'3\n 0\n 1\n'); % number-of-integer-tags, timestep, # of fields to visualize
fprintf(fID,'%d\n',length(hexK)); % < integer-tag > # of elements to view
for ee = 1:length(hexK)
    e = hexK(ee); NpK = NpH;
    fprintf(fID,'%d %d ',e, NpK); %elm-number number-of-nodes-per-element
    for i = 1:NpK
        fprintf(fID,'%f ', pplot(i,e)); %elm-number number-of-nodes-per-element value
    end
    fprintf(fID,'\n');    
end
fprintf(fID,'$EndElementNodeData\n');

% wedge
fprintf(fID,'$ElementNodeData\n');
fprintf(fID,'2\n"wedgeData"\n"Nodal_Poly_Interpolation"\n'); % #-string-tags, "string-tag", interp scheme
fprintf(fID,'1\n 0.1\n'); % number-of-real-tags, < real-tag > - fake time tag
fprintf(fID,'3\n 0\n 1\n'); % number-of-integer-tags, timestep, # of fields to visualize
fprintf(fID,'%d\n',length(wedgK)); % < integer-tag > # of elements to view
for ee = 1:length(wedgK)
    e = wedgK(ee); NpK = NpW;
    fprintf(fID,'%d %d ',e, NpK); %elm-number number-of-nodes-per-element
    for i = 1:NpK
        fprintf(fID,'%f ', pplot(i,e)); %elm-number number-of-nodes-per-element value
    end
    fprintf(fID,'\n');    
end
fprintf(fID,'$EndElementNodeData\n');

% pyr
fprintf(fID,'$ElementNodeData\n');
fprintf(fID,'2\n"pyrData"\n"Nodal_Poly_Interpolation"\n'); % #-string-tags, "string-tag", interp scheme
fprintf(fID,'1\n 0.1\n'); % number-of-real-tags, < real-tag > - fake time tag
fprintf(fID,'3\n 0\n 1\n'); % number-of-integer-tags, timestep, # of fields to visualize

fprintf(fID,'%d\n',length(pyrK)); % < integer-tag > # of elements to view
for ee = 1:length(pyrK)
    e = pyrK(ee); NpK = NpP;
    fprintf(fID,'%d %d ',e, NpK); %elm-number number-of-nodes-per-element
    for i = 1:NpK
        fprintf(fID,'%f ', pplot(i,e)); %elm-number number-of-nodes-per-element value
    end
    fprintf(fID,'\n');    
end
fprintf(fID,'$EndElementNodeData\n');

% tet
fprintf(fID,'$ElementNodeData\n');
fprintf(fID,'2\n"tetData"\n"Nodal_Poly_Interpolation"\n'); % #-string-tags, "string-tag", interp scheme
fprintf(fID,'1\n 0.1\n'); % number-of-real-tags, < real-tag > - fake time tag
fprintf(fID,'3\n 0\n 1\n'); % number-of-integer-tags, timestep, # of fields to visualize
fprintf(fID,'%d\n',length(tetK)); % < integer-tag > # of elements to view
for ee = 1:length(tetK)
    e = tetK(ee); NpK = NpT;
    fprintf(fID,'%d %d ',e, NpK); %elm-number number-of-nodes-per-element
    for i = 1:NpK
        fprintf(fID,'%f ', pplot(i,e)); %elm-number number-of-nodes-per-element value
    end
    fprintf(fID,'\n');
end
fprintf(fID,'$EndElementNodeData\n');

fclose(fID);

disp('outputed p_output.msh')