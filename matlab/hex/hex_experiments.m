% hex experiments
N = 1;
[r s t] = hex_nodes(N);
invV = inv(hex_basis(N,r,s,t));

[rq sq tq w] = hex_cubature(N);
Vq = hex_basis(N,rq,sq,tq);

[rp sp tp] = meshgrid(linspace(-1,1,30));
rp = rp(:); sp = sp(:); tp = tp(:); 
Vp = hex_basis(N,rp,sp,tp);

r1 = [    -1    -1     1     1    -1    -1     1     1]'; 
s1 = [    -1     1    -1     1    -1     1    -1     1]';
t1 = [    -1    -1    -1    -1     1     1     1     1]';
VX = r1 + .5*randn(size(r1)); 
VY = s1 + .5*randn(size(s1)); 
VZ = t1 + .5*randn(size(t1)); 

[xp,yp,zp,~,~,~,~,~,~,~,~,~,Jp] = hex_geom_factors(VX,VY,VZ,rp,sp,tp);
[xq,yq,zq,~,~,~,~,~,~,~,~,~,Jq] = hex_geom_factors(VX,VY,VZ,rq,sq,tq);

Vq = Vq*invV;
M = Vq'*spdiag(w.*Jq(:))*Vq;
%%
color_line3(rp,sp,tp,Jp,'.')
