function [] = testDerivs(N)
% these flags need to be set to "trick" the mass and stiffness computations
ghostBasis = 1;
closure = ghostBasis;
Neumann = 1;

%% discretization parameters
p = 3;
nGhost = (p-1)/2;

%% first set up the overlapping grid
xa = -1;
xb = 1;

%xa = pi/2;
%xb = 3*pi/2;
%N = 80;
x = linspace( xa,xb,N );
dx = x(2)-x(1);

x = [xa-(1:nGhost)*dx,x,xb+(1:nGhost)*dx];
NTot = length(x);

% setup matrices
M = setupMassMatrix( NTot,dx,Neumann,Neumann,p,closure );
K = setupKXXMatrix(  NTot,dx,Neumann,Neumann,p,closure );

if( 1 == 1 )
  k = pi;
  u   = cos(k*x');
  uxx_ex = -k^2*cos(k*x');
elseif( 1 ==  0 )
  u = 1*zeros(size(x'));
  uxx_ex = 0*u;
elseif( 1 == 0 )
  u = ((x+1).^2.*(x-1).^2)';
  uxx_ex = (12*x.^2-4)';
elseif( 1 == 1 )
  u = (x.^3/3-x)';
  uxx_ex = (2*x)';
else
  u = 1*ones(size(x))';
  uxx_ex = dx*ones(size(x))';
end


uxx = M\(K*u);
%uxx = (K*u);
%uxx = M*u;

figure
subplot(2,1,1)
plot( x,uxx,'bx', x,uxx_ex,'k.' );
subplot(2,1,2)
plot( x,uxx-uxx_ex,'k' );

fprintf( '%e\n', max(abs(uxx-uxx_ex)) );

figure
plot( x,M*u,'x' );
%plot( x,u,'b' );
