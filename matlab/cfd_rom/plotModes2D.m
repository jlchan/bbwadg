clear
load Usnap_euler2d_wall.mat

Usnap = reshape(Usnap,size(Usnap,1)/4,size(Usnap,2)*4);
K = sqrt(size(Usnap,1));
[x y] = meshgrid(linspace(-1,1,K));
[V,S,~] = rsvd(Usnap,2*size(Usnap,2));
Usnap = reshape(Usnap,size(Usnap,1)*4,size(Usnap,2)/4);

%%
i = 18;
surf(x,y,reshape(Usnap(1:K^2,i),K,K));shading interp;view(2);axis equal;axis off

%%
% K = 5
e = ones(K-1,1);
Q = diag(e,1)-diag(e,-1);
% Q(1,:) = 0; Q(:,1) = 0;
% Q(end,:) = 0; Q(:,end) = 0;
% Q(1,2) = 1; 
Q(1,end) = -1; 
% Q(end,end-1) = -1; 
Q(end,1) = 1;
% Q(1,2) = -1; Q(K,K) = 1;
Q = .5*Q;


surf(x,y,reshape(V(1:K^2,10),K,K)*Q')
shading interp
view(2)
axis equal
axis off