function layer_driver

clearvars -global

WedgeGlobals

% [VX VY VZ] = wedge_nodes(1); EToV = 1:6; K = 1; EToF = 1:5; EToE = ones(1,5);
% VZ(1) = VZ(1) + .4;
% VZ(6) = VZ(6) + .25;
% VZ = VZ + .2*randn(size(VZ));
% load cube1
load cube2
% load layer_mesh1

N = 3;

WedgeStartUp
% drawWedgeMesh(VX,VY,VZ,EToV);axis equal;view(3)
% keyboard

%% check eigenvalues

if 1
    A = zeros(4*Np*K);
    
    U = zeros(4*Np,K);
    ids = 1:Np*K;   
    
    for i = 1:4*Np*K
        
        U(i) = 1;
        p = reshape(U(ids),Np,K);
        u = reshape(U(ids+Np*K),Np,K);
        v = reshape(U(ids+2*Np*K),Np,K);
        w = reshape(U(ids+3*Np*K),Np,K);
        U(i) = 0;        
        
        [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);
        rhsAll = [rhsp(:); rhsu(:); rhsv(:); rhsw(:)];
        A(:,i) = rhsAll(:);
                
    end
    A(abs(A)<1e-8) = 0;
    [W D] = eig(A); lam = diag(D);
    hold on;plot(lam,'o')
    hold on;plot(1i*linspace(min(imag(lam)),max(imag(lam)),100),'k--','linewidth',2)
    keyboard
else
%     drawWedgeMesh
end


%% solver

uex = @(x,y,z,t) cos(pi/2*x).*cos(pi/2*y).*cos(pi/2*z)*cos(sqrt(3)/2*pi*t);

FinalTime = .1;

% compute time step size
CT = [ 9.926135933275306  18.563576705381948  29.030325215439685  42.988345972840094 58.802145509223905  78.006158337859318  99.271490513773074];
dt = 0.95/max(CT(N)*Fscale(:))

% set initial conditions
p = uex(x,y,z,0); u = zeros(Np,K); v = zeros(Np,K); w = zeros(Np,K); % zero velocities
% [~, bad_id] = max(real(lam));
% W = reshape(W(:,bad_id),4*Np,K); % try unstable mode
% p = W(ids,:); u = W(ids+Np,:); v = W(ids+2*Np,:); w = W(ids+3*Np,:);

% % plot initial cond
% up = Vp*reshape(p,Np,K);
% ids = yp>0;
% h = color_line3(xp(ids),yp(ids),zp(ids),up(ids),'.');
% set(h,'markersize',32);
% keyboard

% outer time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

resp = zeros(Np,K);
resu = zeros(Np,K);
resv = zeros(Np,K);
resw = zeros(Np,K);

% while(tstep<1)
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
    
%     U = [p(:); u(:); v(:); w(:)];
%     U = expm(time*A)*U;
%     U = reshape(U,4*Np,K); 
%     p = U(ids,:); u = U(ids+Np,:); v = U(ids+2*Np,:); w = U(ids+3*Np,:);
    
    % low storage RK
    for INTRK = 1:5        
        [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);

%         U = [p(:);u(:);v(:);w(:)];
%         U(:) = A*U(:);
%         U = reshape(U,Np*K,4);        
%         rhsp = reshape(U(:,1),Np,K); 
%         rhsu = reshape(U(:,2),Np,K); 
%         rhsv = reshape(U(:,3),Np,K); 
%         rhsw = reshape(U(:,4),Np,K); 
        
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
        w = w + rk4b(INTRK)*resw;        
    end;
    
    time = time+dt;
    tstep = tstep+1; % Increment time
    
    if 0 %(mod(tstep,10)==0)
        clf
        up = Vp*reshape(p,Np,K);        
        
        pids = zp < 0;
        h = color_line3(xp(pids),yp(pids),zp(pids),up(pids),'.');
        set(h,'markersize',32);
        title(sprintf('time t = %f',time))
        caxis([-1,1])
        colorbar
        axis([min(xp(:)) max(xp(:)) min(yp(:)) max(yp(:)) min(zp(:)) max(zp(:))])
        view(3);
        
        drawnow
    end
    if mod(tstep,floor(totalSteps/10))==0
        disp(sprintf('on timestep %d / %d\n',tstep,totalSteps));        
    end
    pq = Vq*p;
    pnorm(tstep) = sum(wJ(:).*pq(:).^2);
    
end

% figure
% plot(abs(pnorm))
err = (Vq*p-uex(xq,yq,zq,time)).^2;
err = sqrt(sum(wJ(:).*err(:)));

% hold on
% ids = ones(size(zq(:)));
% h = color_line3(xq(ids),yq(ids),zq(ids),err(ids),'.');
% set(h,'markersize',32);
% title(sprintf('Error at time t = %f = %g',time,err))
% colorbar
% % axis([min(xp(:)) max(xp(:)) min(yp(:)) max(yp(:)) min(zp(:)) max(zp(:))])
% view(3);
err


function [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w)

% Purpose  : Evaluate RHS flux in 3D DG

WedgeGlobals;

Drp = Dr*p; Dsp = Ds*p; Dtp = Dt*p;
Dru = Dr*u; Dsu = Ds*u; Dtu = Dt*u;
Drv = Dr*v; Dsv = Ds*v; Dtv = Dt*v;
Drw = Dr*w; Dsw = Ds*w; Dtw = Dt*w;

dpdx = rx.*Drp + sx.*Dsp + tx.*Dtp;
dpdy = ry.*Drp + sy.*Dsp + ty.*Dtp;
dpdz = rz.*Drp + sz.*Dsp + tz.*Dtp;

dudx = rx.*Dru + sx.*Dsu + tx.*Dtu;
dvdy = ry.*Drv + sy.*Dsv + ty.*Dtv;
dwdz = rz.*Drw + sz.*Dsw + tz.*Dtw;
divU = dudx + dvdy + dwdz;

p_jump = zeros(Nfp,K);
nu_jump = zeros(Nfp,K);

p_jump(:) = p(vmapP)-p(vmapM);
p_jump(mapB) = -2*p(vmapB); % boundary conditions

u_jump  = u(vmapP) - u(vmapM);
v_jump  = v(vmapP) - v(vmapM);
w_jump  = w(vmapP) - w(vmapM);
nu_jump(:) = nx(:).*u_jump + ny(:).*v_jump + nz(:).*w_jump;
nu_jump(mapB) = 0; % extra neumann BCs??

tau_p = 1;
tau_u = 1;
flux_p = .5*(p_jump*tau_p - nu_jump);%.*Fscale(:); %Fscale built into stored LIFTs?
flux_u = .5*(nu_jump*tau_u - p_jump);%.*Fscale(:); %

rhsp = -divU;
rhsu = -dpdx;
rhsv = -dpdy;
rhsw = -dpdz;

flux_p = reshape(flux_p,Nfp,K);
flux_u = reshape(flux_u,Nfp,K);
for e = 1:K
    rhsp(:,e) = rhsp(:,e) + LIFT{e}*(flux_p(:,e));
    rhsu(:,e) = rhsu(:,e) + LIFT{e}*(nx(:,e).*flux_u(:,e));
    rhsv(:,e) = rhsv(:,e) + LIFT{e}*(ny(:,e).*flux_u(:,e));
    rhsw(:,e) = rhsw(:,e) + LIFT{e}*(nz(:,e).*flux_u(:,e));
end


