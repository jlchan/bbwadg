function Elasticity2D

% clear all, clear
clear -global *

Globals2D

K1D = 16;
N = 7;
c_flag = 0;
FinalTime = 1;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
VX = (VX+1)/2; VY = (VY+1)/2;
StartUp2D;

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%%
global Nfld mu lambda Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = ones(size(x));
lambda = ones(size(x));

% mu = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% lambda = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% ids = mean(y) < .5;
% mu(:,ids) = 2*mu(:,ids);
% lambda(:,ids) = 2*lambda(:,ids);
% mu = Vq*mu;
% lambda = Vq*lambda;

tau = [ones(2,1); ones(3,1)/4];
% tau = ones(Nfld,1);

useWADG = 0;
%%

if 0
    tauvec = 1;
    for ii = 1:length(tauvec)
        tau = tauvec(ii)*[ones(2,1); ones(3,1)/4];
%         tau = tauvec(ii)*[1 1 1 1 1];
        u = zeros(Nfld*Np*K,1);
        rhs = zeros(Nfld*Np*K,1);
        A = zeros(Nfld*Np*K);
        ids = 1:Np*K;
        for i = 1:Nfld*Np*K
            u(i) = 1;
            for fld = 1:Nfld
                U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
            end
            rU = ElasRHS2D(U);
            u(i) = 0;
            for fld = 1:Nfld
                rhs(ids + (fld-1)*Np*K) = rU{fld};
            end
            A(:,i) = rhs(:);
            if (mod(i,100)==0)
                disp(sprintf('on col %d out of %d\n',i,Np*K*Nfld))
            end
        end
        lam = eig(A);        
        points{ii} = [real(lam) imag(lam)];
%         plot(lam,'.','markersize',32)
%         hold on
%         title(sprintf('Largest real part = %g\n',max(real(lam))))
%         axis equal
%         drawnow
        max(abs(lam))
    end
    keyboard

    %% track eigs
    max_linking_distance = 10;
    max_gap_closing = Inf;
    debug = true;
    disp('running simple tracker')
    [ tracks, adjacency_tracks ] = simpletracker(points,...
        'MaxLinkingDistance', max_linking_distance, ...
        'MaxGapClosing', max_gap_closing, ...
        'Debug', debug);
    
    figure
    plot(points{1}(:,1),points{1}(:,2),'bo','markersize',10,'linewidth',2,'DisplayName','\tau = 0')
    hold on
    [~,id] = min(abs(tauvec-1/2));
    plot(points{id}(:,1),points{id}(:,2),'r^','markersize',10,'DisplayName','\tau = 1.0')
    plot(points{end}(:,1),points{end}(:,2),'gs','markersize',10,'DisplayName',sprintf('\\tau = %1.1f',tauvec(end)))
    
    % % exact eigs
    % lam_ex = pi/2*(1:Np*K); lam_ex = [lam_ex -lam_ex];
    % plot(zeros(size(lam_ex)),lam_ex,'s','markersize',11,'DisplayName','Exact eigenvalues')
    
    legend show
    n_tracks = numel(tracks);
    all_points = vertcat(points{:});
    for i_track = 1 : n_tracks
        
        % We use the adjacency tracks to retrieve the points coordinates. It
        % saves us a loop.
        
        track = adjacency_tracks{i_track};
        track_points = all_points(track, :);
        
        plot(track_points(:,1), track_points(:, 2), 'k','linewidth',2)
    end
    axis equal
    set(gca,'fontsize',15)
    grid on
    %print(gcf,'-dpng','trackedEigs.png')
    keyboard
end

%% params setup
k = 1; % frequency of solution
W = (2*k-1)/2*pi;
x0 = mean(VX); y0 = mean(VY) + .1;
p = exp(-50^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);

w = sqrt(2);
u1 = @(x,y,t) cos(w*pi*t)*cos(pi*x).*sin(pi*y);
u2 = @(x,y,t) -cos(w*pi*t)*sin(pi*x).*cos(pi*y);

if 0
    for t = 0:.1:1;
        clf
        vv = Vp*u1(x,y,t);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis([0 1 0 1 -1 1])
        drawnow
    end
    return
end

U{1} = u1(x,y,0);
U{2} = u2(x,y,0);
U{3} = u;
U{4} = u;
U{5} = u;

U{1} = u;
U{2} = p;
U{3} = u;
U{4} = u;
U{5} = u;

%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 1/(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
figure
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        rhs = ElasRHS2D(U);
        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:Nfld
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;
    
    if mod(tstep,10)==0
        clf
        
        %         subplot(1,2,1)
        p = U{2}; % trace(S)
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        %         axis([0 1 0 1 -1 1])
        %
        %         subplot(1,2,2)
        %         p = U{2};
        %         vv = Vp*p;
        %         color_line3(xp,yp,vv,vv,'.');
        %         axis tight
        
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

p = U{1}; % trace(S)
vv = Vp*p;
err = (abs(u1(xp,yp,FinalTime)-vv));
figure
color_line3(xp,yp,err,err,'.')



function [rhs] = ElasRHS2D(U)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda Vq Pq tau useWADG

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u;
    us = Ds*u;
    Ux{fld} = rx.*ur + sx.*us;
    Uy{fld} = ry.*ur + sy.*us;
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% Impose reflective boundary conditions
%dU{1}(mapB) = -2*U{1}(vmapB);
%dU{2}(mapB) = -2*U{2}(vmapB);
% dU{5}(mapB) = -2*U{5}(vmapB)*0;

% evaluate upwind fluxes
fc{1} = nx.*dU{3} + ny.*dU{5};
fc{2} = nx.*dU{5} + ny.*dU{4};
fc{3} = dU{1}.*nx;
fc{4} = dU{2}.*ny;
fc{5} = dU{2}.*nx + dU{1}.*ny;

% penalization terms
fp{1} = dU{1} + nx.*ny.*dU{2};
fp{2} = dU{2} + nx.*ny.*dU{1};
fp{3} = nx.*(nx.*dU{3} + ny.*dU{5});
fp{4} = ny.*(nx.*dU{5} + ny.*dU{4});
fp{5} = dU{5} + nx.*ny.*(dU{3} + dU{4});

% zero traction BCs - set nx*Sxx + ny*Sxy = nx*Sxy + ny*Syy = 0
fc{1}(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
fc{2}(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
fp{3}(mapB) = -2*nx(mapB).*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
fp{4}(mapB) = -2*ny(mapB).*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));

for fld = 1:Nfld
    flux{fld} = tau(fld)*fp{fld} - fc{fld};
end

% compute right hand sides of the PDE's
rr{1} =  -divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  -divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  -du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  -du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  -du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu/2]

if useWADG
    for fld = 3:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = Pq*((mu/2) .* rr{5});
else
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu/2) .* rr{5};
end

return;

