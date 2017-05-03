% function circle_inclusion
clear

rad = .25;
N = 16;
h = 2/N;

theta = linspace(-1+2/N, 1, N)';
cnodes = rad*[cos(pi*theta),sin(pi*theta)];

% plot(cnodes(:,1),cnodes(:,2),'o');  axis equal;return

%  Here are all the nodes/edges.
%
node = [
    -1.0, -1.0;
    1.0, -1.0;
    1.0, 1.0;
    -1.0, 1.0];

edge = [
    1,2;
    2,3;
    3,4;
    4,1];

nc = size(cnodes,1);
cedges = zeros(nc,2);
cedges(:,1) = 1:nc;
cedges(:,2) = [2:nc 1];
cedges = length(node) + cedges;

node = [node;cnodes];
edge = [edge;cedges];

%  Face 1 is the outer strip.
%  Face 2 is the inner square.

face{1} = 1:length(node);
face{2} = 4 + (1:length(cnodes));

%  HDATA allows us to set options.
%    HDATA.HMAX is the maximum triangle side length.

hdata = [];
hdata.hmax = h;
% hdata.edgeh = [face{2}(:) ones(size(face{2}(:)))];
%
%  MESHFACES output is:
%
%    P     = Nx2 array of nodal XY co-ordinates.
%    T     = Mx3 array of triangles as indicies into P, defined with a
%            counter-clockwise node ordering.
%    FNUM  = Mx1 array of face numbers for each triangle in T.
%
[ p, t, fnum, stats ] = meshfaces ( node, edge, face, hdata );

% fnum==1 and fnum==2 give different faces
circelems = find(fnum==2);
for ee = 1:length(circelems)
    e = circelems(ee);
    tri = t(e,:);
    for i = 1:3
        pt = p(t(e,i),:);
        
        for j = 1:length(cnodes)
            if j == 1
                id1 = length(cnodes);
                id2 = j;                
            else
                id1 = j-1;
                id2 = j;
            end
            a = norm(cnodes(id1,:)-pt);
            b = norm(cnodes(id2,:)-pt);
            c = norm(cnodes(id1,:)-cnodes(id2,:));
            tol = 1e-6;
            if a < tol || b < tol
                % do nothing
            elseif  abs(a+b-c) < tol
                pt = rad*(pt)/norm(pt);
                p(t(e,i),:) = pt;                                
            end
        end
    end
end

%  [p,t1] = smoothmesh(p,t(fnum==1,:),1,0.1);
%  t(fnum==1,:) = t1;
%  [p,t] = fixmesh(p,t);
% [p,t2] = smoothmesh(p,t(fnum==2,:),10,0.1);
% t(fnum==2,:) = t2;

clf
% figure('Name','Mesh')
plot(p(:,1),p(:,2),'b.','markersize',1)
hold on;
% Colour mesh for each face
col = ['b','r','g','o','m'];
for k = 1:length(face)
    colk = mod(k,length(col));
    if (colk==0)
        colk = length(col);
    end
    patch('faces',t(fnum==k,:),'vertices',p,'facecolor','w','edgecolor',col(colk));
end
patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')

axis on
