function M = setupMassMatrix( Ntot,h,iBCLeft,iBCRight,p,closure )
M = zeros( Ntot,Ntot );

ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;
difference    = 4;

Neumann = 1;
Dirichlet = 2;

% determine gauss quadrature to integrate exactly 2p order polynomial
[qx,qw] = lgwt( p+1,0,1 );
phii = 0*qx;
phij = 0*qx;

nGhost = (p-1)/2;
ileft  = 1+nGhost;
iright = Ntot-nGhost;

if( 1 == 1 )
    %% perform quadratures for a single canonical element and store the results
    %ints = zeros( 2*p,2*p );
    %i = 0;
    %for j = i-p:i+p
    %  %% integrate phi_i against phi_j
    %  for k = -(p+1)/2:(p+1)/2-1
    %    xa = k;
    %    xb = k+1;
    %    for iq = 1:length(phii);
    %      phii(iq) = phi( xa+qx(iq),p );
    %      phij(iq) = phi( xa+(i-j)+qx(iq),p );
    %    end
    %    ints( j-(i-p)+1,k+(p+1)/2+1 ) = sum( (phii.*phij).*qw );
    %  end
    %end
    fname = sprintf( 'integrationTables/massMatrix_p%i.dat', p );
    fid = fopen( fname,'r' );
    ints = reshape( fread( fid,'real*8' ),2*p+1,p+1 );
    fclose( fid );
    
    %% now use these results to form mass matrix
    for i = 1:Ntot
        %% ith equation
        for j = max(i-p,1):min(i+p,Ntot)
            %% integrate phi_i against phi_j
            tmp = 0;
            for k = -(p+1)/2:(p+1)/2-1
                if( i+k >= ileft & i+k+1 <= iright )
                    tmp = tmp+ints(j-(i-p)+1,k+(p+1)/2+1);
                end
            end
            M(i,j) = h*tmp;
        end
    end
    
else
    
    
    for i = 1:Ntot
        %% ith equation
        for j = max(i-p,1):min(i+p,Ntot)
            %% integrate phi_i against phi_j
            %fprintf( '%i %i \n\n', i,j );
            tmp = 0;
            for k = -(p+1)/2:(p+1)/2-1
                xa = k;
                xb = k+1;
                if( i+k >= ileft & i+k+1 <= iright )
                    %fprintf( 'integrating [%i %i]\n',xa+i,xb+i );
                    for iq = 1:length(phii);
                        phii(iq) = phi( xa+qx(iq),p );
                        phij(iq) = phi( xa+(i-j)+qx(iq),p );
                    end
                    tmp = tmp+sum( (phii.*phij).*qw );
                end
            end
            %      pause
            M(i,j) = h*tmp;
        end
    end
end


if( closure == compatibility )
    if( iBCLeft == Neumann )
        for ig = 1:nGhost
            M(:,ileft+ig)  = M(:,ileft+ig) +M(:,ileft-ig);
            M(ileft+ig,:)  = M(ileft+ig,:) +M(ileft-ig,:);
        end;
    else % Dirichlet
        for ig = 1:nGhost
            M(:,ileft+ig)  = M(:,ileft+ig) -M(:,ileft-ig);
            M(ileft+ig,:)  = M(ileft+ig,:) -M(ileft-ig,:);
        end;
        ileft = ileft+1;
    end;
    
    if( iBCRight == Neumann )
        for ig = 1:nGhost
            M(:,iright-ig) = M(:,iright-ig)+M(:,iright+ig);
            M(iright-ig,:) = M(iright-ig,:)+M(iright+ig,:);
        end;
    else % Dirichlet
        for ig = 1:nGhost
            M(:,iright-ig) = M(:,iright-ig)-M(:,iright+ig);
            M(iright-ig,:) = M(iright-ig,:)-M(iright+ig,:);
        end;
        iright = iright-1;
    end
    M = M(ileft:iright,ileft:iright);
    
elseif( closure == ghostBasis )
    inds = [1:ileft-1];
    if( iBCLeft ~= Dirichlet )
        inds = [inds,ileft];
    end
    inds = [inds,ileft+1:iright-1];
    if( iBCRight ~= Dirichlet )
        inds = [inds,iright];
    end
    inds = [inds,iright+1:Ntot];
    M = M(inds,inds);
    
elseif( closure == extrapolation )
    nCompat = 0;
    for ig = nGhost:-1:1
        ec = getExtrapCoeffs( p+1 );
        if( ig <= nCompat & iBCLeft == Neumann )
            %if( ig == nGhost & iBCLeft == Neumann )
            ec = zeros( 2*ig+1,1 );
            ec(end) = 1;
        elseif( ig <= nCompat & iBCLeft == Dirichlet )
            %elseif( ig ==nGhost & iBCLeft == Dirichlet )
            ec = zeros( 2*ig+1,1 );
            ec(end) = -1;
        end
         % ileft  = 1+nGhost;
        iGhost = ileft-ig;
        for k = 2:length(ec)
            M(:,iGhost+k-1) = M(:,iGhost+k-1)+ec(k)*M(:,iGhost);
            M(iGhost+k-1,:) = M(iGhost+k-1,:)+ec(k)*M(iGhost,:);
        end
    end
    
    for ig = nGhost:-1:1
        ec = getExtrapCoeffs( p+1 );
        if( ig <= nCompat & iBCRight == Neumann )
            %if( ig == nGhost & iBCRight == Neumann )
            ec = zeros( 2*ig+1,1 );
            ec(end) = 1;
        elseif( ig <= nCompat & iBCRight == Dirichlet )
            %elseif( ig == nGhost & iBCRight == Dirichlet )
            ec = zeros( 2*ig+1,1 );
            ec(end) = -1;
        end
        iGhost = iright+ig;
        for k = 2:length(ec)
            M(:,iGhost-k+1) = M(:,iGhost-k+1)+ec(k)*M(:,iGhost);
            M(iGhost-k+1,:) = M(iGhost-k+1,:)+ec(k)*M(iGhost,:);
        end
    end
    
    if( iBCLeft == Dirichlet )
        ileft = ileft+1;
    end
    
    if( iBCRight == Dirichlet )
        iright = iright-1;
    end
    M = M(ileft:iright,ileft:iright);
    
else % closure == difference
    for ig = nGhost:-1:1
        if( iBCLeft == Neumann )
            dc = generalDifferenceWeights( 0,[-ig:-ig+p+1],p+1 );
            %dc = generalDifferenceWeights( 0,[-ig:-ig+p+1+ig],p+1+ig );
            ec = dc(:,2);
            ec = -ec/(ec(1));
            %pause
            %ec = zeros( 2*ig+1,1 );
            %ec(end) = 1;
        elseif( iBCLeft == Dirichlet )
            ec = getExtrapCoeffs( p+1 );
            %ec = zeros( 2*ig+1,1 );
            %ec(end) = -1;
        end
        iGhost = ileft-ig;
        for k = 2:length(ec)
            M(:,iGhost+k-1) = M(:,iGhost+k-1)+ec(k)*M(:,iGhost);
            M(iGhost+k-1,:) = M(iGhost+k-1,:)+ec(k)*M(iGhost,:);
        end
    end
    
    for ig = nGhost:-1:1
        if( iBCRight == Neumann )
            dc = generalDifferenceWeights( 0,[-ig:-ig+p+1],p+1 );
            %dc = generalDifferenceWeights( 0,[-ig:-ig+p+1+ig],p+1+ig );
            ec = dc(:,2);
            ec = -ec/(ec(1));
            %ec = zeros( 2*ig+1,1 );
            %ec(end) = 1;
        elseif( iBCRight == Dirichlet )
            %         dc = generalDifferenceWeights( 0,[-ig:-1,1:-ig+p+1],p );
            %         ec = dc(:,1);
            %         ec = [ec(1:ig);1;ec(ig+1:end)]
            %         ec = -ec/(ec(1))
            ec = getExtrapCoeffs( p+1 );
            %         pause
            %         ec = zeros( 2*ig+1,1 );
            %         ec(end) = -1;
        end
        iGhost = iright+ig;
        for k = 2:length(ec)
            M(:,iGhost-k+1) = M(:,iGhost-k+1)+ec(k)*M(:,iGhost);
            M(iGhost-k+1,:) = M(iGhost-k+1,:)+ec(k)*M(iGhost,:);
        end
    end
    
    if( iBCLeft == Dirichlet )
        ileft = ileft+1;
    end;
    
    if( iBCRight == Dirichlet )
        iright = iright-1;
    end
    M = M(ileft:iright,ileft:iright);
    %%%%
    
end