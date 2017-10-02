clear

%% compute reference solution

FinalTime = .1;

a = 1;
cfun = @(x,y) 1 + .5*sin(a*pi*x).*sin(a*pi*y); % smooth velocity
% cfun = @(x,y) 1 + .5*sin(2*pi*x + 1).*sin(pi*y + 2);
[x y] = meshgrid(linspace(-1,1,100));
% pcolor(x,y,cfun(x,y));shading interp; colorbar; return

Nref = 50;
[pref invVquad] = spectral_wave(Nref,cfun,FinalTime); % N = 100 for reference
[x y] = meshgrid(JacobiGL(0,0,Nref));

% surf(x,y,pref);shading interp;return

% %% check quadrature dependence
% ref = 4;
% K1D = 2.^ref;
% N = 4;
% [rq sq w] = Cubature2D(3*N);
% [r s] = Nodes2D(N); [r s] = xytors(r,s);
% Vq = Vandermonde2D(N,rq,sq)/Vandermonde2D(N,r,s);
% [p0 x y] = save_solution_driver(N,K1D,cfun,FinalTime,0);
% [p1 x y] = save_solution_driver(N,K1D,cfun,FinalTime,1);
% 
% xq = Vq*x; yq = Vq*y;
% 
% Vqref = Vandermonde2DQuad(Nref,xq(:),yq(:))*invVquad;
% prefq = reshape(Vqref*pref(:),size(xq,1),size(xq,2));
% 
% J = 2*.5.^(2*ref); % for uniform mesh
% errWDG = sqrt(sum(sum(J*spdiag(w)*(prefq - Vq*p0).^2)));
% errDG = sqrt(sum(sum(J*spdiag(w)*(prefq - Vq*p1).^2)));
% format shorte
% errDG 
% errWDG
% format
% return
% %% compute

%%
N = 4;
[rq sq w] = Cubature2D(2*N+1);
[r s] = Nodes2D(N); [r s] = xytors(r,s);
Vq = Vandermonde2D(N,rq,sq)/Vandermonde2D(N,r,s);
    
sk = 1;
for ref = 1:4
    disp(sprintf('on ref %d\n',ref))
    K1D = 2.^ref;
    [p0 x y] = save_solution_driver(N,K1D,cfun,FinalTime,0);
%     [p1 x y] = save_solution_driver(N,K1D,cfun,FinalTime,1);
    
    xq = Vq*x; yq = Vq*y;
    
    Vqref = Vandermonde2DQuad(Nref,xq(:),yq(:))*invVquad;        
    prefq = reshape(Vqref*pref(:),size(xq,1),size(xq,2));
    
    J = 2*.5.^(2*ref); % for uniform mesh
    err0(sk) = sqrt(sum(sum(J*spdiag(w)*(prefq - Vq*p0).^2)));
%     err1(sk) = sqrt(sum(sum(J*spdiag(w)*(prefq - Vq*p1).^2)));
    sk = sk + 1;    
end
'hi'
% plot3(xq,yq,prefq)
% figure;plot3(xq,yq,Vq*p0)
%%

h = 2*.5.^(1:length(err0)); 
loglog(h,err0,'o-');hold on
loglog(h,1e-2*h.^(N+1),'--')
return

%%

h = 2*.5.^(1:length(err0)); 
h = h(:); err0 = err0(:); err1 = err1(:);

print_pgf_coordinates(h,err0)
print_pgf_coordinates(h,err1)
format shorte
fprintf('%4.4e & ',err0);fprintf('\n')
fprintf('%4.4e & ',err1);fprintf('\n')
format

ids = ref-1:ref;
C1 = [h(ids).^0 log(h(ids))]\log(err0(ids));
C2 = [h(ids).^0 log(h(ids))]\log(err1(ids));

rates = [C1(2), C2(2)]
return

%%

% err0 = [   0.011496577446261   0.002109218124822   0.000117841156514 0.000005758942629];
% err1 = [   0.011226029605943   0.002097372904872   0.000118290263828   0.000005795756681];
loglog(h,err0,'o-');hold on
loglog(h,err1,'s-')

ids = ref-1:ref;
C0 = [h(ids).^0 log(h(ids))]\log(err0(ids));
C1 = [h(ids).^0 log(h(ids))]\log(err1(ids));
title(sprintf('Rates = %f and %f, expected rate %f',C0(2),C1(2),N+1))
legend('Low storage projection','c-dependent mass matrix')
set(gca,'fontsize',14)