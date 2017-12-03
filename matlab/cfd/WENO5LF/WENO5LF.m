%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 1-D Euler system of equations with 5th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%        dq_i/dt + df_i/dx = 0, for x \in [a,b] and i =1,. ..,D
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2012.08.25
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Sod's shock tube problem (IC=1)
%
% t=0                                 t=tEnd
% Density                             Density
%   ****************|                 *********\
%                   |                           \
%                   |                            \
%                   |                             ****|
%                   |                                 |
%                   |                                 ****|
%                   ***************                       ***********
%
% coded by Manuel A. Diaz, 2012.12.27. Last modif: 29.04.2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: C.-W. Shu, High order weighted essentially non-oscillatory schemes
% for convection dominated problems, SIAM Review, 51:82-126, (2009). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 
% 1. A fully conservative finite volume implementation of the method of
% lines (MOL) using WENO5 associated with SSP-RK33 time integration method. 
% 2. Sharpenning of contact discontinuities is NOT implemented here.

clear; %close all; clc;
global gamma

%% Parameters
CFL     = 0.5;	% CFL number
tFinal	= 1.8;	% Final time
nE      = 10001;  % Number of cells/Elements
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = 01;	% 10 IC cases are available
plot_fig= 1;

% Discretize spatial domain
a=-5; b=5; dx=(b-a)/nE; nx=nE+1; x=linspace(a,b,nx);

% Set IC
[rho0,u0,p0] = Euler_IC1d(x,IC);

opt=2;
if opt==1
    % sod initial condition
    x0 = (b+a)/2;
    rho0 = (1*(x < x0) + .125*(x >= x0));
    p0 = (1*(x < x0) + .1*(x >= x0));
    u0 = 0*x;
elseif opt==2
    rhoRex = @(x) 1 + .2*sin(5*x);
    rhoL = 3.857143; rhoR = rhoRex(x(end));
    uL   = 2.629369; uR = 0;
    pL   = 10.3333;  pR = 1;
    
    rhoRex = @(x) 2 + .2*sin(5*x);
    rhoL = 4; rhoR = rhoRex(x(end));
    uL   = 2; uR = 0;
    pL   = 10; pR = 2;
        
    rho0 = rhoL*(x < -4) + (rhoRex(x)).*(x >= -4);
    u0 = uL*(x < -4);
    p0 = pL*(x < -4) + pR*(x >= -4);
end
% plot(x,p0,'o')
% return
% E = p/(gamma-1) + .5*rho.*u.^2;
% m = rho.*u;

% plot(x,rho0,'o')
% return
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % Total Energy density
a0 = sqrt(gamma*p0./rho0);            % Speed of sound
q0=[rho0; rho0.*u0; rho0.*E0];        % vec. of conserved properties

% Exact solution
[xe,rhoe,ue,pe,ee,te,Me,se] = ...
   EulerExact(rho0(1),u0(1),p0(1),rho0(nx),u0(nx),p0(nx),tFinal,gamma);

% Discretize time domain
lambda0=max(abs(u0)+a0); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

%% Solver Loop
while t<tFinal
 
    % RK Initial step
    qo = q;
    
    % 1st stage
    dF=WENO5LF1d(lambda,q,dx);     q = qo-dt*dF; 
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs
    
    % 2nd Stage
    dF=WENO5LF1d(lambda,q,dx);     q = 0.75*qo+0.25*(q-dt*dF);
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs

    % 3rd stage
    dF=WENO5LF1d(lambda,q,dx);     q = (qo+2*(q-dt*dF))/3;
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs
   
    % compute primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:)./rho; p=(gamma-1)*rho.*(E-0.5*u.^2);
    a=sqrt(gamma*p./rho); if min(p)<0; error('negative pressure found!'); end
    
    % Update dt and time
    lambda=max(abs(u)+a); dt=CFL*dx/lambda; if t+dt>tFinal; dt=tFinal-t; end
    
    % Update time and iteration counter
	t=t+dt; it=it+1;
    
    % Plot figure
    if rem(it,10) == 0
        if plot_fig == 1;
            subplot(2,2,1); plot(x,rho,'.b'); legend('\rho')
            subplot(2,2,2); plot(x,u,'.m'); legend('u')
            subplot(2,2,3); plot(x,p,'.k'); legend('p')
            subplot(2,2,4); plot(x,E,'.r'); legend('E')
            title(sprintf('time = %f',t))
            
        end
	drawnow
    end
end

% Calculation of flow parameters
a = sqrt(gamma*p./rho); M = u./a; % Mach number [-]
p_ref = 101325;             % Reference air pressure (N/m^2)
rho_ref= 1.225;             % Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho)); 
                            % Entropy w.r.t reference condition
ss = log(p./rho.^gamma);    % Dimensionless Entropy
Q = rho.*u;                 % Mass Flow rate per unit area
e = p./((gamma-1)*rho);     % internal Energy

return
%% Final plot
offset=0.05;
s1=subplot(2,3,1); plot(x,rho,'or',xe,rhoe,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); plot(x,u,'or',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
s3=subplot(2,3,3); plot(x,p,'or',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); plot(x,ss,'or',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); plot(x,M,'or',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); plot(x,e,'or',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'MOL WENO-LF Euler Solver');