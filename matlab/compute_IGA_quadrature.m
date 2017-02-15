% compute IGA quad rules
function compute_IGA_quadrature

global reflect f NT tt DBr t NB Vq wq

NB = 4;
Ksub = NB;
N = NB+Ksub-1;
[Nv, VX, Ksub, EToV] = MeshGen1D(-1,1,Ksub);
map = @(r) reshape((1/Ksub)*(repmat(r,1,Ksub)+1) + ...
    repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);

% local quadrature
[rq, wq] = JacobiGQ(0,0,2*NB+2); % overkill gauss quadrature
rq = map(rq);
wq = repmat(wq,1,Ksub)*(1/Ksub);
wq = wq(:);

t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
Vq = bspline_basismatrix(NB+1,t,rq);
M = Vq'*diag(wq)*Vq;
M(abs(M)<1e-8) = 0;

% plot(rq,rq*0,'o')
% hold on
% plot(VX,VX*0,'s','markersize',18)
% plot(rq,Vq)
% return

%% compute quad rules

NT = 2*NB; % target order - for GQ, 2*NB, for GLL 2*NB-1
tt = sort([VX(1)*ones(1,NT+1) repmat(VX(2:end-1),1,NB+1) VX(end)*ones(1,NT+1)]); % product of splines
% tt = sort([t t]);
Bq = bspline_basismatrix(NT+1,tt,rq);

if 0
    plot(rq,Bq)
    for i = 1:size(Vq,2)
        for j = 1:size(Vq,2)
            vv = Vq(:,i).*Vq(:,j);
            % check representation in basis
            qerr = norm(Bq*(Bq \ vv) - vv);
            if qerr > 1e-10
                qerr
                keyboard
            end
            
        end
    end
    
    disp(sprintf('quad pts + weights = %d, size of product space = %d, piecewise GQ = %d\n',...
        (N+1),ceil(size(Bq,2)/2),(NB+1)*Ksub)) % dim of space, # qpts to integrate target spline space, # gq pts
end

f = Bq'*wq; % exact integrals
Nq = ceil(size(Bq,2)/2)-1;

if 0    
    rprev = spline_quadrature(NB-1);
    reprev = linspace(-1,1,length(rprev));
    
    % bootstrap up
    rq0 = (Vandermonde1D(N,linspace(-1,1,Nq+1))/Vandermonde1D(N,reprev))*rprev;
    
    if mod(length(rq0),2)==0 % if #pts even, reflect
        reflect = @(x) sort([-x(:); x(:)]); % reflects points, map [0,1] to [-1,1]
        %     reflect = @(x) sort([-1;-x(:); x(:); 1]); % include endpoints for GLL
    else % if #pts = odd, assume 0 is the middle point
        reflect = @(x) sort([-x(:); 0; x(:)]); % reflects points, map [0,1] to [-1,1]
        %     reflect = @(x) sort([-1;-x(:);0; x(:); 1]); % include endpoints for GLL
    end
    
    tol = 1e-10;
    opt = optimoptions('fsolve','FunctionTolerance',tol,'StepTolerance',tol,'OptimalityTolerance',tol); %,'SpecifyObjectiveGradient',true);
    rq0reflect = rq0((ceil(length(rq0)/2)+1):end);
    r = reflect(fsolve(@(x) residual(reflect(x)),rq0reflect,opt));
end

% just use predefined points
for i = 1:N+1
    r(i) = mean(t((i+1):(i+NB))); % greville
end
% r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

% recover weights
Vq = bspline_basismatrix(NT+1,tt,r);
w = (Vq') \ f;

norm(residual(r))

Vq = bspline_basismatrix(NB+1,t,rq);
V = bspline_basismatrix(NB+1,t,r);
M = Vq'*diag(wq)*Vq;
Mh = V'*diag(w)*V;
norm(M-Mh,'fro')

lam = eig(Mh,M);
[min(lam),max(lam)]

disp(sprintf('Nqpts = %d, size of space = %d\n',Nq+1,size(Bq,2)))

keyboard

rp = linspace(-1,1,250)';
Vp = bspline_basismatrix(NB+1,t,rp);
W = inv(V);
plot(rp,Vp*W);
hold on
plot(r,r*0,'ks')


function [F J] = residual(rq) % r01 = half of points on (0,1)

global f NT tt DBr

Bq = bspline_basismatrix(NT+1,tt,rq);

% solve for weights
wq = (Bq') \ f;

% residual
F = Bq'*wq - f;

norm(F)
