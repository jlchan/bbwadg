function [uerrOut,verrOut,h,tRun] = wave( N,p,closure,plotOption,nPeriods,tf )

%% p sets the reconstruction order
%% iclosure sets the type of boundary closure
%% nPeriods is approximately the number of periods in the domain

ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;
difference    = 4;

warning( 'off','MATLAB:nearlySingularMatrix' )

icOption = 6;
%icOption = 3;

Neumann = 1;
Dirichlet = 2;
iBCLeft  = 0;
iBCRight = 0;

%nPeriods = 12;
%nPeriods = 97;
%nPeriods = 194;

iBCLeft  = Neumann;
iBCRight = Neumann;
if( icOption == 6 )
    iBCLeft  = Neumann;
    iBCRight = Neumann;
elseif( icOption == 7 )
    iBCLeft  = Neumann;
    iBCRight = Dirichlet;
elseif( icOption == 8 )
    iBCLeft  = Dirichlet;
    iBCRight = Neumann;
elseif( icOption == 9 )
    iBCLeft  = Dirichlet;
    iBCRight = Dirichlet;
end

t = 0;
%tf = 1.25;
c = 1;
if( closure == ghostBasis | extrapolation | closure == difference )
    cfl = 0.25;
else
    cfl = 0.9;
end

uem = [];
vem = [];
tem = [];


if( icOption == 6 )
    xa = -nPeriods+1/2;
    xb = nPeriods-1/2;
elseif( icOption == 7 )
    xa = -nPeriods+1/2;
    xb = nPeriods;
elseif( icOption == 8 )
    xa = -nPeriods;
    xb = nPeriods-1/2;
else
    xa = -nPeriods;
    xb = nPeriods;
end

h = (xb-xa)/(N-1);
nGhost = (p-1)/2;
Ntot = N+2*nGhost;
x = linspace( xa-nGhost*h,xb+nGhost*h,Ntot );
u = zeros(2,Ntot);
ileft  = 1+nGhost;
iright = Ntot-nGhost;

if( icOption <= 4 | icOption >= 6 )
    for i = 1:Ntot
        ul = icfunc( x(i)-c*t,icOption );
        ur = icfunc( x(i)+c*t,icOption );
        u(1,i) = 0.5*(ul+ur);
    end
else
    for i = 1:Ntot
        u(1,i) = sin(pi*(x(i)-c*t));
        u(2,i) = -pi*c*cos(pi*(x(i)-c*t));
    end
end

M   = setupMassMatrix( Ntot,h,iBCLeft,iBCRight,p,closure );
Kxx = setupKXXMatrix(  Ntot,h,iBCLeft,iBCRight,p,closure );

%cond(M)
%pause

dt_temp = cfl*h/abs(c)*2.78/pi;
nt = ceil( tf/dt_temp );
dt = tf/nt;
%fprintf( 'time steps: %i\n', nt );

%% here we optionally lump the mass matrix
if( 1 == 1 )
    [nm,~] = size(M);
    MDiag  = diag(sum(M'));
    pskip = p+2;
    inds = pskip:nm-(pskip)+1;
    M(inds,:) = MDiag(inds,:);
    %     M(:,inds) = MDiag(:,inds);
end

[L,U] = lu(sparse(M),0);

tic
for step = 1:nt
    %% 2*p+2 order
    %% this is a bit dangerous since I do not know if the time integrator includes
    %%   the imaginary axis for arbitrary p
    u1 = u;
    up0 = u;
    for ns = 1:2*p+2
        %for ns = 1:p+1
        %for ns = 1:4
        up = wave_step( up0,M,Kxx,c,iBCLeft,iBCRight,p,L,U,closure );
        u = u+dt^(ns)/factorial(ns)*up;
        up0= up;
    end
    
    t = t+dt;
    
    if( plotOption ~= 0 && (mod(step,25)==0 || step==nt))
        uex = 0*x;
        vex = 0*x;
        if( icOption <= 4 | icOption >= 6 )
            for i = 1:Ntot
                ul = icfunc( x(i)-c*t,icOption );
                ur = icfunc( x(i)+c*t,icOption );
                uex(i) = 0.5*(ul+ur);
                
                vl = icfuncp( x(i)-c*t,icOption );
                vr = icfuncp( x(i)+c*t,icOption );
                vex(i) = 0.5*(c*vr-c*vl);
            end
        else
            for i = 1:Not
                uex(i) = sin(pi*(x(i)-c*t));
                vex(i) = -pi*c*cos(pi*(x(i)-c*t));
            end
        end
        
        if( plotOption == 1 )
            %figure(1)
            %             subplot(2,1,1), plot( x(ileft:iright),u(1,ileft:iright),'bx-', x(ileft:iright),uex(ileft:iright),'k-' );
            %             subplot(2,1,2), plot( x(ileft:iright),u(2,ileft:iright),'bx-', x(ileft:iright),vex(ileft:iright),'k-' );
            subplot(2,1,1), plot( x(ileft:iright),u(1,ileft:iright),'bx-');
            subplot(2,1,2), plot( x(ileft:iright),u(2,ileft:iright),'bx-');
        elseif( plotOption == 2 )
            %figure(1)
            %ileft = 1;
            %iright = length(x);
            subplot(4,1,1), plot( x(ileft:iright),u(1,ileft:iright),'bx-', x(ileft:iright),uex(ileft:iright),'k-' );
            subplot(4,1,2), plot( x(ileft:iright),u(2,ileft:iright),'bx-', x(ileft:iright),vex(ileft:iright),'k-' );
            subplot(4,1,3), plot( x(ileft:iright),u(1,ileft:iright)-uex(ileft:iright),'kx-' );
            subplot(4,1,4), plot( x(ileft:iright),u(2,ileft:iright)-vex(ileft:iright),'kx-' );
        elseif( plotOption == 3 )
            tem = [tem,t];
            uem = [uem,max(abs(u(1,ileft:iright)-uex(ileft:iright)))];
            vem = [vem,max(abs(u(2,ileft:iright)-vex(ileft:iright)))];
            plot( tem,uem,'bx', tem,vem,'rs' );
        end
        title(sprintf('time = %f\n',step*dt))
        drawnow
%           pause(.01);
        %pause
    end
    
end

tRun = toc;
%%%%%
%%%%%

uex = 0*x;
vex = 0*x;
if( icOption <= 4 | icOption >= 6 )
    for i = 1:Ntot
        ul = icfunc( x(i)-c*t,icOption );
        ur = icfunc( x(i)+c*t,icOption );
        uex(i) = 0.5*(ul+ur);
        
        vl = icfuncp( x(i)-c*t,icOption );
        vr = icfuncp( x(i)+c*t,icOption );
        vex(i) = 0.5*(c*vr-c*vl);
    end
else
    for i = 1:Ntot
        uex(i) = sin(pi*(x(i)-c*t));
        vex(i) = -pi*c*cos(pi*(x(i)-c*t));
    end
end

uerr = abs(uex(ileft:iright)-u(1,ileft:iright));
verr = abs(vex(ileft:iright)-u(2,ileft:iright));

uerrL1 = sum(uerr)/N;
verrL1 = sum(verr)/N;

uerrL2 = sqrt(sum(uerr.^2)/N);
verrL2 = sqrt(sum(verr.^2)/N);

uerrMax = max(uerr);
verrMax = max(verr);

uerrOut = [uerrL1;uerrL2;uerrMax];
verrOut = [verrL1;verrL2;verrMax];

%fprintf( '%e %e %e\n',   sum(uerr)/N,sqrt(sum(uerr.^2)/N),max(uerr) );
%fprintf( '%e %e %e\n\n', sum(verr)/N,sqrt(sum(verr.^2)/N),max(verr) );

%if( plotOption == 0 )
%  figure
%  subplot(2,1,1), plot( x,u(1,:)-uex,'bx-' );
%  subplot(2,1,2), plot( x,u(2,:)-vex,'bx-' );
%end
