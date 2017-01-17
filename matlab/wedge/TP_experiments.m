clear
u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
VX = v; VY = w; VZ = u; % flipping coordinates for Gmsh
a = .5;
VX = VX + a*randn(size(VX)); VY = VY + a*randn(size(VX)); VZ = VZ + a*randn(size(VX));
VX = VX(:); VY = VY(:); VZ = VZ(:);

N = 1;
[r s t w] = wedge_cubature(N);
[x,y,z,~,~,~,~,~,~,~,~,~,J] = wedge_geom_factors(VX,VY,VZ,r,s,t);

A(:,1) = ones(size(r));
A(:,2) = r;
A(:,3) = s;
A(:,4) = t;
A(:,5) = r.*s;
A(:,6) = s.*t;

C = A\J;
norm(J-A*C)
%%

NcTri = length(r)/(N+1);
J = reshape(J,NcTri,N+1);
r = reshape(r,NcTri,N+1);
t = reshape(t,NcTri,N+1);
% plot3(r,s,t,'.')
% text(r,s,t,num2str((1:length(r))'))

%%
w = w.*J;


sk = 1;
for j = 0:N
    for i = 0:N        
        for k = 0:N-i
            ids(sk,:) = [j i k sk];
            sk = sk + 1;
        end
    end
end
ids

% regular operators
[V Vr Vs Vt] = wedge_basis(N,r,s,t);

% TP operators
[V2 V2r V2t V1 V1s NpTri NcTri] = wedge_TP_basis(N);

% NpTri = (N+1)*(N+2)/2;
% NcTri = size(V2,1);
u = randn(size(V,2),1);
uu = reshape(u,NpTri,N+1);

% interp to cub
uc = zeros(NcTri,N+1);
for i = 1:N+1
    uc(:,i) = (V2*uu(:,i));
end
for j = 1:NcTri
    uc(j,:) = V1*(uc(j,:)');
end
norm(V*u - uc(:))

% check one entry
Vu = V*u;
uc = V1*(uu');
uc = V2*(uc);

return
% r deriv
dudr = zeros(NcTri,N+1);
for i = 1:N+1
    dudr(:,i) = (V2r*uu(:,i));
end
for j = 1:NcTri
    dudr(j,:) = dudr(j,:)*V1';
end
norm(Vr*u - dudr(:))

% s deriv
duds = zeros(NcTri,N+1);
for i = 1:N+1
    duds(:,i) = (V2*uu(:,i));
end
for j = 1:NcTri
    duds(j,:) = duds(j,:)*V1s';
end
norm(Vs*u - duds(:))

% t deriv
dudt = zeros(NcTri,N+1);
for i = 1:N+1
    dudt(:,i) = (V2t*uu(:,i));
end
for j = 1:NcTri
    dudt(j,:) = V1*(dudt(j,:)');
end
norm(Vt*u - dudt(:))

w = reshape(w,NcTri,N+1);
dudr = w.*dudr;
VUr = zeros(NpTri,N+1);
for i = 1:N+1
    VUr(:,i) = V2'*dudr(:,i);
end
for j = 1:NpTri
    VUr(j,:) = V1'*(VUr(j,:)');
end
norm(V'*(w(:).*(Vr*u))-VUr(:))
