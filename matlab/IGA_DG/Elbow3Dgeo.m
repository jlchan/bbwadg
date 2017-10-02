function [x y z geofacs] = Elbow3Dgeo(r,s,t)
% function Elbow3Dgeo(r,s,t)

% h = 1; w = 1; rad = .5; R = 1; T = rad/2;
h = 2; w = 2; rad = 0.5; R = 2; T = 0.2; 

[N,N_2,N_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3,P,W] = NURBS_Pipe_Elbow_Volume(h,w,rad,R,T);

d = 3;
[K,C_operators,IEN] = Extract_Basis(N,N_2,N_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3);

[P_b,w_b] = Extract_Geometry(d,K,C_operators,IEN,P,W);

% Bezier_Plotter(N,N_2,N_3,K,P_b,w_b);
% return

if nargin==0
    %r1D = linspace(-1,1,N+1)';
    r1D = linspace(-1,1,8)';
    [r s t] = meshgrid(r1D); r = r(:); s = s(:); t = t(:);
end

Np = (N+1)^3;
cx = reshape(P_b(:,1,:),Np,K);
cy = reshape(P_b(:,2,:),Np,K);
cz = reshape(P_b(:,3,:),Np,K);
cw = reshape(w_b,Np,K);
% cw = cw.^0; % restore polynomial approx
% keyboard


if 0
    [V1D Vr1D] = bern_basis_1D(N,r1D);
    V = kron(kron(V1D,V1D),V1D);
    Vr = kron(kron(V1D,Vr1D),V1D);
    Vs = kron(kron(V1D,V1D),Vr1D);
    Vt = kron(kron(Vr1D,V1D),V1D);
else
    [V1Da Vr1Da] = bern_basis_1D(N,r);
    [V1Db Vr1Db] = bern_basis_1D(N,s);
    [V1Dc Vr1Dc] = bern_basis_1D(N,t);
    V = zeros(length(r),(N+1)^3);
    Vr = zeros(length(r),(N+1)^3);
    Vs = zeros(length(r),(N+1)^3);
    Vt = zeros(length(r),(N+1)^3);
    sk = 1;
    for k = 1:N+1
        for j = 1:N+1
            for i = 1:N+1
                V(:,sk) = V1Da(:,i).*V1Db(:,j).*V1Dc(:,k);
                Vr(:,sk) = Vr1Da(:,i).*V1Db(:,j).*V1Dc(:,k);
                Vs(:,sk) = V1Da(:,i).*Vr1Db(:,j).*V1Dc(:,k);
                Vt(:,sk) = V1Da(:,i).*V1Db(:,j).*Vr1Dc(:,k);
                sk = sk + 1;
            end            
        end
    end    
    
end


w = V*cw;
x = (V*(cx.*cw))./w; 
y = (V*(cy.*cw))./w;
z = (V*(cz.*cw))./w;
wr = Vr*cw;
ws = Vs*cw;
wt = Vt*cw;
xr = (Vr*(cx.*cw))./w - wr.*x./w; 
xs = (Vs*(cx.*cw))./w - ws.*x./w;
xt = (Vt*(cx.*cw))./w - wt.*x./w;

yr = (Vr*(cy.*cw))./w - wr.*y./w; 
ys = (Vs*(cy.*cw))./w - ws.*y./w;
yt = (Vt*(cy.*cw))./w - wt.*y./w;

zr = (Vr*(cz.*cw))./w - wr.*z./w; 
zs = (Vs*(cz.*cw))./w - ws.*z./w;
zt = (Vt*(cz.*cw))./w - wt.*z./w;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
% rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
% sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
% tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

rxJ =  (ys.*zt - zs.*yt); ryJ = -(xs.*zt - zs.*xt); rzJ = (xs.*yt - ys.*xt);
sxJ = -(yr.*zt - zr.*yt); syJ =  (xr.*zt - zr.*xt); szJ = -(xr.*yt - yr.*xt);
txJ =  (yr.*zs - zr.*ys); tyJ = -(xr.*zs - zr.*xs); tzJ = (xr.*ys - yr.*xs);

geofacs.rxJ = rxJ; geofacs.sxJ = sxJ; geofacs.txJ = txJ;
geofacs.ryJ = ryJ; geofacs.syJ = syJ; geofacs.tyJ = tyJ;
geofacs.rzJ = rzJ; geofacs.szJ = szJ; geofacs.tzJ = tzJ;
geofacs.J = J;


