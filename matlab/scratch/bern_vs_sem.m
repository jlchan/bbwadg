
N = 4;
[r1D w1D] = JacobiGQ(0,0,N); 
% [r1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,r1D); w1D = sum(inv(V*V'),2);

[VB VrB] = bern_basis_1D(N,r1D); D = VB\VrB;

[r s t] = meshgrid(r1D); r = r(:); s = s(:); t = t(:); 
[wr ws wt] = meshgrid(w1D); w = wr(:).*ws(:).*wt(:);
[V Vr Vs Vt] = bern_hex(N,r,s,t);


[rf sf] = meshgrid(r1D); rf = rf(:); sf = sf(:); tf = -ones(size(rf));
[wrf wsf] = meshgrid(w1D); wf = wrf(:).*wsf(:);

% compute lift matrix
V = bern_hex(N,r,s,t);
Vf = bern_hex(N,rf,sf,tf); 

% nodal basis
[rL sL tL] = meshgrid(JacobiGL(0,0,N)); % GLL nodes
rL = rL(:); sL = sL(:); tL = tL(:);

VN = nodal_hex(N,rL,sL,tL);
V = nodal_hex(N,r,s,t)/VN;
Vf = nodal_hex(N,rf,sf,tf)/VN;

% [rp sp tp] = meshgrid(linspace(-1,1,15)); rp = rp(:); sp = sp(:); tp = tp(:);
% Vp = nodal_hex(N,rp,sp,tp)/VN;
% for i = 1:(N+1)^3
%     h = color_line3(rp,sp,tp,Vp(:,i),'.');set(h,'markersize',32);
%     hold on;plot3(rf,sf,tf,'ro','markersize',16)    
%     hold on;plot3(r,s,t,'go','markersize',18)
%     view(3)
%     pause
% end
%%

M = V'*diag(w)*V;
   
Mf = Vf'*diag(wf)*Vf;

LIFT = M\Mf;
LIFT(abs(LIFT)<1e-6) = 0; 
% LIFTf = LIFT(:,1:(N+1)^2);

imagesc(LIFT)

figure
% plot(VX([1 2 4 3 1]),VY([1 2 4 3 1]),'o-');hold on
plot3(x,y,z,'o');hold on
plot3(xf,yf,zf,'r*')
quiver3(xf,yf,zf,nx,ny,nz)
% axis equal


