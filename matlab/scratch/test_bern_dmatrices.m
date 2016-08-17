N = 2;
[a b c] = meshgrid(linspace(-1,1-(1/(N+1)),N+1));
a = a(:); b = b(:); c = c(:);
r = .5*(1+a).*(1-c)-1;
s = .5*(1+b).*(1-c)-1;
t = c;

% plot3(r,s,t,'.')
[V Vr Vs Vt Va Vb Vc] = bern_pyr(N,r,s,t);
u = V\s;

norm(s-V*u,'fro')
norm(1-Vs*u,'fro')

Dr = V\Vr;
Dr(abs(Dr)<1e-8) = 0;


for j = 1:size(V,2)
    Sa(:,j) = (2./(1-c)).*V(:,j);
end

% for j = 1:size(V,2)
%     clf
%     color_line3(r,s,t,V(:,j),'.')
%     view(3)
%     pause
% end