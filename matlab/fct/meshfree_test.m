clear
N = 6;
load ~/Downloads/QuadratureRules.mat

quadpts = Q_GaussLegendre{N};
% quadpts = Q_GaussLobatto{N};
xy = quadpts.Points;
w = quadpts.Weights;

x = xy(:,1);
y = xy(:,2);

% face nodes
[r1D w1D] = JacobiGQ(0,0,N);
ef = ones(size(r1D));
xf = [r1D; r1D; -ef];
yf = [-ef; -r1D; r1D];
wf = [w1D; sqrt(2)*w1D; w1D];
nx = [0*ef; ef; -ef];
ny = [-ef; ef; 0*ef];

fid = [];
for i = 1:length(xf)
    fid = [fid; find(abs(x-xf(i))+abs(y-yf(i))<1e-10)];
end

ep = 3/(N+1);
[xe ye] = rstoxy(x,y);

tc = linspace(0,2*pi,100)';
xc = cos(tc); yc = sin(tc);

% plot(x,y,'o')
% hold on
% for i = 1:10:length(x(:))    
%     sc = ep;
%     plot(x(i)+sc*xc,y(i)+sc*yc,'k-')    
% end

adj = zeros(length(x(:)));
for i = 1:length(x(:))
    d2 = (xe(i)-xe).^2 + (ye(i) - ye).^2;
    adj(i,d2 < ep*ep) = 1;
    adj(i,i) = 0;
end

% % compute adjacency using equilateral tri
% plot(xe,ye,'bo','linewidth',2,'markersize',16)
% hold on
% for i = 1:length(x(:))
%     id = find(adj(i,:));    
%     for j = 1:length(id)        
%         plot([xe(i) xe(id(j))],[ye(i) ye(id(j))],'k-')
%     end   
% end

% basis for p1
Vfun = @(x,y) [x(:).^0 x(:).^1 y(:).^1];
dxVfun = @(x,y) [0*x(:) x(:).^0 0*x(:)];
dyVfun = @(x,y) [0*x(:) 0*x(:) x(:).^0];

Lx = {};
for k = 1:3
    Lx{k} = zeros(length(x(:)));
    f{k} = 0*x(:);
end
for i = 1:length(x(:))
    id = find(adj(i,:));
    for j = 1:length(id)
        xa = .5*(x(i) + x(id(j)));
        ya = .5*(y(i) + y(id(j)));
        phia = Vfun(xa,ya);         
        for k = 1:3
            Lx{k}(i,i) = Lx{k}(i,i) - phia(k);
            Lx{k}(i,id(j)) = phia(k);
        end
    end
    
    dphi = dxVfun(x(i),y(i));    
    for k = 1:3
        f{k}(i) = f{k}(i) + w(i)*dphi(k);
        if ismember(i,fid)            
            idf = find(i==fid);
            phif = Vfun(x(i),y(i));
            f{k}(i) = f{k}(i) + wf(idf)*phif(k)*nx(idf); 
        end
    end
    
end

% for k = 1:3
%     mu{k} = pinv(Lx{k})*f{k};
% end
% 
% 
% for i = 1:length(x(:))    
%     id = find(adj(i,:));
%     for j = 1:length(id)
%         xa = .5*(x(i) + x(id(j)));
%         ya = .5*(y(i) + y(id(j)));
%         phia = Vfun(xa,ya);
%         for k = 1:3
%             muscaled{k}(i,j) = (mu{k}(j)-mu{k}(i))*phia(k)*Cjk;
%         end        
%     end
% end

% reconstruct
u = sin(.5*x+y);

% build virtual points
xa = [];
ya = [];
for i = 1:length(x)
    idi = find(adj(i,:));
    xai = .5*(x(i) + x(idi)); 
    yai = .5*(y(i) + y(idi));
    xa = [xa; xai];
    ya = [ya; yai];
end

ur = 0*xa;
offset = 0;
for i = 1:length(x)
    idi = find(adj(i,:));
    id = [i idi];
    Vi = Vfun(x(id),y(id));
    idf = offset + (1:length(idi)); 
    ur(idf) = ur(idf) + Vfun(xa(idf),ya(idf)) * (Vi \ u(id));    
    offset = offset + length(idi);    
end

clf
vv = u;
color_line3(x,y,vv,vv,'.')
hold on
vv = ur;
color_line3(xa,ya,vv,vv,'o')
    