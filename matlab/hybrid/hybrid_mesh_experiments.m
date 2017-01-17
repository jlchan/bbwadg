% script to check face numbering of each element
N = 3;

%%
hex_mesh
K = size(EToV,1);
minJ = zeros(K,1);
for e = 1:K
    v = EToV(e,:); v = v(v > 0);
    minJ(e) = min_hex_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
end
disp(sprintf('min jacobian = %f\n',min(minJ)))

e = 1;
v = EToV(e,:); v = v(v > 0);

hold on
plot3(VX(v),VY(v),VZ(v),'.')
text(VX(v)+.05,VY(v)+.05,VZ(v)+.05,num2str((1:8)'))
text(mean(VX(v)),mean(VY(v)),mean(VZ(v)),sprintf('ELEM %d',e))

ids = [1 2 3 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 1')
ids = [1 2 6 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 2')
ids = [2 3 7 6]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 3')
ids = [3 4 8 7]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 4')
ids = [4 1 5 8]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 5')
ids = [5 6 7 8]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 6')

% fidsQ{1} = [1 4 3 2];
% fidsQ{2} = [1 2 6 5];
% fidsQ{3} = [2 3 7 6];
% fidsQ{4} = [3 4 8 7];
% fidsQ{5} = [4 1 5 8];
% fidsQ{6} = [5 6 7 8];
% v = EToV(e,:); v = v(v > 0);
% [rf sf tf wf] = surface_cubature(N,VX(v),VY,VZ,fidsQ);
% color_line3(rf,sf,tf,wf,'.')

%%
prism_mesh2;
K = size(EToV,1);
minJ = zeros(K,1);
for e = 1:K
    v = EToV(e,:); v = v(v > 0);
    minJ(e) = min_wedge_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
end
disp(sprintf('min jacobian = %f\n',min(minJ)))

e = 1;
v = EToV(e,:); v = v(v > 0);

plot3(VX(v),VY(v),VZ(v),'.')
text(VX(v)+.05,VY(v)+.05,VZ(v)+.05,num2str((1:6)'))
text(mean(VX(v)),mean(VY(v)),mean(VZ(v)),sprintf('ELEM %d',e))

hold on

ids = [1 2 3]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 1')
ids = [4 5 6]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 2')
ids = [1 2 5 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 4')
ids = [3 1 4 6]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 3')
ids = [2 3 6 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 5')

fidsW{1} = [1 2 3];
fidsW{2} = [4 5 6];
fidsW{3} = [3 1 4 6];
fidsW{4} = [1 2 5 4];
fidsW{5} = [2 3 6 5];

% v = EToV(e,:); v = v(v > 0);
% 
% [rf sf tf wf] = surface_cubature(N,VX,VY,VZ,v,fidsW);
% 
% color_line3(rf,sf,tf,wf,'.')

% min_wedge_jacobian(VX(v),VY(v),VZ(v));
%%
pyr_mesh;
K = size(EToV,1);
minJ = zeros(K,1);
for e = 1:K
    v = EToV(e,:); v = v(v > 0);
    minJ(e) = min_pyr_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
end
disp(sprintf('min jacobian = %f\n',min(minJ)))

e = 1;
v = EToV(e,:); v = v(v > 0);
% VZ(v) = -VZ(v);
hold on
plot3(VX(v),VY(v),VZ(v),'.')
text(VX(v)+.05,VY(v)+.05,VZ(v)+.05,num2str((1:5)'))
text(mean(VX(v)),mean(VY(v)),mean(VZ(v)),sprintf('ELEM %d',e))

ids = [1 2 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 1')
ids = [2 3 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 3')
ids = [3 4 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 4')
ids = [4 1 5]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 2')
ids = [1 2 3 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 5')

fidsP{1} = [1 2 5];
fidsP{2} = [2 3 5];
fidsP{3} = [3 4 5];
fidsP{4} = [4 1 5];
fidsP{5} = [1 2 3 4];

v = EToV(e,:); v = v(v > 0);
[rf sf tf wf] = surface_cubature(N,VX,VY,VZ,v,fidsP);

color_line3(rf,sf,tf,wf,'.')

% min_pyr_jacobian(VX(v),VY(v),VZ(v))

%%
tet_mesh;

K = size(EToV,1);
minJ = zeros(K,1);
for e = 1:K
    v = EToV(e,:); v = v(v > 0);
    minJ(e) = min_tet_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
end
disp(sprintf('min jacobian = %f\n',min(minJ)))
figure
% get vertex nodal bases
e = 1;
v = EToV(e,:); v = v(v > 0);

hold on
plot3(VX(v),VY(v),VZ(v),'.')
text(VX(v)+.05,VY(v)+.05,VZ(v)+.05,num2str((1:4)'))
text(mean(VX(v)),mean(VY(v)),mean(VZ(v)),sprintf('ELEM %d',e))

ids = [1 2 3]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 1')

ids = [1 2 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 2')

ids = [2 3 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 4')
%
ids = [3 1 4]; pv = v([ids ids(1)]);
plot3(VX(pv),VY(pv),VZ(pv),'linewidth',2)
text(mean(VX(pv)),mean(VY(pv)),mean(VZ(pv)),'face 3')

% min_tet_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
fidsT{1} = [1 2 3];
fidsT{2} = [1 2 4];
fidsT{3} = [2 3 4];
fidsT{4} = [3 1 4];

v = EToV(e,:); v = v(v > 0);
[rf sf tf wf] = surface_cubature(N,VX,VY,VZ,v,fidsT);

color_line3(rf,sf,tf,wf,'.')

%% hybrid mesh

hybrid_mesh
K = size(EToV,1);
minJ = zeros(K,1);
type = zeros(K,1);
for e = 1:K
    v = EToV(e,:); v = v(v > 0); type(e) = nnz(v);
    if e==1 || e==2
        plot3(VX(v),VY(v),VZ(v),'.');hold on
        text(VX(v)+.05,VY(v)+.05,VZ(v)+.05,num2str((1:type(e))'))
        text(mean(VX(v)),mean(VY(v)),mean(VZ(v)),sprintf('ELEM %d',e))
    end

    switch type(e)
        case 4
            minJ(e) = min_tet_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
        case 5
            minJ(e) = min_pyr_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
        case 6 
            minJ(e) = min_wedge_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
        case 8
            minJ(e) = min_hex_jacobian(VX(v),VY(v),VZ(v)); % to check ordering
    end    
        
end
disp(sprintf('min jacobian = %f\n',min(minJ)))

%% compute surface Jacobians


