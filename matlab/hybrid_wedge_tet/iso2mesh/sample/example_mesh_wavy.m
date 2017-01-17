x = linspace(-1,1,50);
y = linspace(-1,1,100);
z = linspace(-1,1,25);
[x y z] = meshgrid(x,y,z);

% generate random wavy layers
Nw = 4;
c1 = randn(Nw)+.25; c2 = randn(Nw)+.25;
wavy1 = zeros(size(x)); wavy2 = zeros(size(x));
for ii = 1:Nw
    for jj = 1:Nw
        wavy1 = wavy1+c1(ii,jj)*sin(ii*x+.1).*sin(jj*y+.3)/(ii*jj);
        wavy2 = wavy2+c2(ii,jj)*sin(ii*x+.2).*sin(jj*y+.4)/(ii*jj);
    end
end
layer1 = .25*wavy1;
sec1 = z < layer1;
%vel = (sec1 + sec2 + 1);
vel = (sec1 + 1);
% vel(:) = 1;
vel=uint8(vel);


% call cgalmesher to mesh the segmented volume
% this will take 30 seconds on an Intel P4 2.4GHz PC

fprintf(1,'meshing layers...\n');

%maxvol = '1=8:2=4:3=2';
maxvol = 8;%'1=16:2=8:3=4';
% maxvol = 50;
opt = 5;
[node,elem,face]=v2m(vel,uniquetol(double(vel(:))),opt,maxvol,'cgalmesh');
  
% node_in = node; face_in = face; elem_in = elem;
% [no2,fc2]=meshcheckrepair(node_in,elem_in);
% [node,face]=meshresample(no2,fc2,0.2);

% %fprintf(1,'remeshing surface...\n');
% [newno,newfc]=remeshsurf(node,face,1);
% newno=sms(newno,newfc(:,1:3),3,0.5);

% laplacian smoothing
nsmooth = 0;
% method = 'laplacian';
method = 'laplacianhc';
% method = 'lowpass';
conn=meshconn(face(:,1:3),size(node,1));
n1=node;
for i= 1:nsmooth
  n1=smoothsurf(n1,[],conn,1,0.5,method);  
end
node = n1;
quality=meshquality(node(:,1:3),elem(:,1:4)); 

figure
% hs=plotmesh(node,face,'x>y');
% hs=plotmesh(node(:,1:3),face);
tag=elem(find(elem(:,5)==2),:);
wmsurf=volface(tag(:,1:4));
hs=plotmesh(node(:,1:3),wmsurf);
axis equal;
% shading interp
title('cross-cut view of the generated surface mesh');

return


% trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3));
return
%%
qe = find(quality<.1);
% qe = qe(1:5);

plotmesh(node(:,1:3),elem(qe(id),:))

for i = 1:length(qe)
    J = -1;
    ids = 1:4;     
    while J<0        
        e = elem(qe(i),ids);
        x = node(e,1); y = node(e,2); z = node(e,3);
        [rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = geofacs(x,y,z);
        J = J(1);                
        ids([2 3]) = ids([3 2]);
    end
    Je(i) = J;    
end

%%
figure
tag=elem(find(elem(:,5)==2),:);
wmsurf=volface(tag(:,1:4));

x = linspace(-1,1,100);
y = linspace(-1,1,200);
[x2D y2D] = meshgrid(x,y);

wavy1 = zeros(size(x2D)); wavy2 = zeros(size(x2D));
for ii = 1:Nw
    for jj = 1:Nw
        wavy1 = wavy1+c1(ii,jj)*sin(ii*x2D+.1).*sin(jj*y2D+.3)/(ii*jj);
        wavy2 = wavy2+c2(ii,jj)*sin(ii*x2D+.2).*sin(jj*y2D+.4)/(ii*jj);
    end
end
layer1 = .2*wavy1 - .3;
layer2 = .2*wavy2 +.3;
mesh(25*layer1'+60)
hold on
mesh(25*layer2'-10)
% return

hs=plotmesh(node(:,1:3),wmsurf);
axis equal;
% shading interp
title('cross-cut view of the generated surface mesh');
