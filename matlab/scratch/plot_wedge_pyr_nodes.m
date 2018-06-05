clear
N = 5;

tol = 1e-7;

if 0
    [r s t] = pyr_nodes(N);
    [r1 s1 t1] = pyr_nodes(1);
    
    fids = cell(5,1);
    fids{1} = abs(r+1)<tol;
    fids{2} = abs(s+1)<tol;
    fids{3} = abs(t+1)<tol;
    fids{4} = abs(s+t)<tol;
    fids{5} = abs(r+t)<tol;    
    
else
    [r s t] = wedge_nodes(N);
    [r1 s1 t1] = wedge_nodes(1);
    % wedge
    fids{1} = abs(s+1)<tol;
    fids{2} = abs(s-1)<tol;
    fids{3} = abs(r+1)<tol;
    fids{4} = abs(t+1)<tol;
    fids{5} = abs(r+t)<tol;
end

% compute edge ids
eids = zeros(size(r));
for f1 = 1:5
    for f2 = 1:5
        if f1~=f2
            eids = (eids | (fids{f1} & fids{f2}));
        end
    end
end
fidsAll = zeros(size(r));
for f = 1:5
    fidsAll = fidsAll | fids{f};
end
fids = fidsAll;
vids = ismembertol(r,r1) & ismembertol(s,s1) &ismembertol(t,t1);

eids = eids ~= vids;
fids = (fids ~= eids) ~= vids;
iids = (((1:length(r)) ~= fids)~=eids)~=vids;

flag =  zeros(size(r));
flag(vids) = 1;
flag(eids) = 4;
flag(fids) = 8;
flag(iids) = 16;

[xs,ys,zs] = sphere(75);
cs = 0*xs;
ra= 0.05;
hold on
for n=1:length(r)
    ha = surf(ra*xs+r(n), ra*ys+s(n), ra*zs + t(n), flag(n)+1 + cs);
    
    shading interp
    material shiny
end

% col(1,:) = [.9 .9 .9];
% col(2,:) = [.8 .8 .8];
% col(3,:) = [.7 .7 .7];
% col(4,:) = [.6 .6 .6];
% col(5,:) = [.5 .5 .5];
% col(6,:) = [.4 .4 .4];
% col(7,:) = [0  0   0];
% colormap(col)

lighting gouraud
camlight

hold off; axis off; axis equal





