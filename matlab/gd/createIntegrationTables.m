for p = 1:2:129

  fprintf( 'working on p=%i\n',p );

  [qx,qw] = lgwt( p+1,0,1 );
  phii = 0*qx;
  phij = 0*qx;

  nGhost = (p-1)/2;

  ints = zeros( 2*p+1,p+1 );

  %%%%%
  %% first mass matrix
  %%%%%
  % determine gauss quadrature to integrate exactly 2p order polynomial

  %% perform quadratures for a canonical element and store the results
  i = 0;
  for j = i-p:i+p
    %% integrate phi_i against phi_j
    for k = -(p+1)/2:(p+1)/2-1
      xa = k;
      xb = k+1;
      for iq = 1:length(phii);
        phii(iq) = phi( xa+qx(iq),p );
        phij(iq) = phi( xa+(i-j)+qx(iq),p );
      end
      ints( j-(i-p)+1,k+(p+1)/2+1 ) = sum( (phii.*phij).*qw );
    end
  end
  fname = sprintf( 'integrationTables/massMatrix_p%i.dat', p );
  fid = fopen( fname,'w+' );
  fwrite( fid,ints,'real*8' );
  fclose( fid );

  %%%%%
  %% now stiffness matrix for u_xx
  %%%%%
  %% perform quadratures for a single canonical element and store the results
  i = 0;
  for j = i-p:i+p
    %% integrate (phi_i)_x against (phi_j)_x
    for k = -(p+1)/2:(p+1)/2-1
      xa = k;
      xb = k+1;
      for iq = 1:length(phii);
        phii(iq) = dphi( xa+qx(iq),p );
        phij(iq) = dphi( xa+(i-j)+qx(iq),p );
      end
      ints( j-(i-p)+1,k+(p+1)/2+1 ) = sum( (phii.*phij).*qw );
    end
  end
  fname = sprintf( 'integrationTables/KXXMatrix_p%i.dat', p );
  fid = fopen( fname,'w+' );
  fwrite( fid,ints,'real*8' );
  fclose( fid );

end

%  fid = fopen( fname,'r' );
%  %A = fread( fid,'real*8' );
%  A = reshape( fread( fid,'real*8' ),2*p+1,p+1 );
%  fclose( fid );
