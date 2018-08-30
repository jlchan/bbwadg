lineWidth = 2;
fs = 14;
ms = 14;

ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;
difference    = 4;

%nPeriods = 50;
nPeriods = 10;
plotOption = 0;
tf = 1.25;

closure = difference;

writePlots = 0;

%res     = 50:1:600;
res     = 100:10:250;
%res     = 100:50:400;
nSims   = length(res);
pp      = 1:2:15;
nOrders = length(pp)
nNorms  = 3;

uerr   = zeros(nSims,nNorms,nOrders);
verr   = zeros(nSims,nNorms,nOrders);
dx     = zeros(nSims,1);
t      = zeros(nSims,nOrders);

for ires = 1:nSims
  Nx = res(ires);
  fprintf( 'sims using %i points\n',Nx );
  ue   = [];
  ve   = [];
  tc   = [];

  for ip = 1:nOrders
    p = pp(ip);
   [ueloc,veloc,dxloc,tRun] = wave( Nx,p,closure,plotOption,nPeriods,tf );
    tc = [tc,tRun];
    ue = [ue,ueloc];
    ve = [ve,veloc];
  end

  t( ires,: )   = tc;

  uerr(ires,:,:)   = ue;
  verr(ires,:,:)   = ve;

  dx(ires) = dxloc;

end


%%%%%
%%%%%
%%%%%         uuuuuuu
%%%%%
%%%%%

%%%%%
%% max norm
%%%%%
dxRef = linspace(min(dx)/1.1,max(dx)*1.1,201);
inorm = 3;
figure
set(gca,'FontSize',fs);
loglog( dx,uerr(:,inorm,1),'bx','lineWidth',lineWidth,'MarkerSize',ms  );
hold on
loglog( dx,uerr(:,inorm,2),'ms','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,3),'go','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,4),'r^','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,5),'ks','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,6),'bd','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,7),'mo','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dx,uerr(:,inorm,8),'g^','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dx,uerr(:,inorm,9),'rs','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dx,uerr(:,inorm,10),'kd','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dx,uerr(:,inorm,11),'bo','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,1e0*dxRef.^2,'b-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,4e0*dxRef.^5,'m-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,1e1*dxRef.^7,'g-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,4e1*dxRef.^8,'r-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,5e2*dxRef.^10,'k-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,5e3*dxRef.^12,'b-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,2e4*dxRef.^14,'m-','lineWidth',lineWidth,'MarkerSize',ms  );
loglog( dxRef,2e5*dxRef.^16,'g-','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dxRef,5e5*dxRef.^18,'r-','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dxRef,8e6*dxRef.^20,'k-','lineWidth',lineWidth,'MarkerSize',ms  );
%loglog( dxRef,8e7*dxRef.^22,'b-','lineWidth',lineWidth,'MarkerSize',ms  );
hold off
axis( [min(dx)/1.1,max(dx)*1.1,1e-16,10] );
xlabel( 'h' );
ylabel( 'error' );
title( 'max norm error in u' );
if( writePlots == 1 )
  plotName = sprintf('images/wave_u_ND_difference_Lm.eps');
  fprintf('Saving file=[%s]\n',plotName);
  print('-depsc2',plotName);
end

