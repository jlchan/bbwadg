function up = wave_step( uold,M,Kxx,c,iBCLeft,iBCRight,p,L,U,closure )

ghostBasis    = 1;
compatibility = 2;
extrapolation = 3;
difference    = 4;

Neumann = 1;
Dirichlet = 2;

[junk,N] = size(uold);
up = 0*uold;

nGhost = (p-1)/2;
ileft  = 1+nGhost;
iright = N-nGhost;

if( closure ~= ghostBasis )
    if( iBCLeft == Dirichlet )
        ileft = ileft+1;
    end
    if( iBCRight == Dirichlet )
        iright = iright-1;
    end
    inds = [ileft:iright];
    
else % compatibility == ghostBasis
    inds = [1:ileft-1];
    if( iBCLeft == Neumann )
        inds = [inds,ileft];
    end
    inds = [inds,ileft+1:iright-1];
    if( iBCRight == Neumann )
        inds = [inds,iright];
    end
    inds = [inds,iright+1:N];
end

%up(1,inds) = uold(2,inds);
up(1,:) = uold(2,:);

if( 1 == 0 )
    b = uold(1,inds)';
    tmp_rhs = 0*b;
    NN = length(b);
    for j = 1:NN
        tmp = 0;
        for k = max(1,j-p):min(NN,j+p)
            tmp = tmp+Kxx(j,k)*b(k);
        end
        tmp_rhs(j) = tmp;
    end
    up(2,inds) = c^2*(U\(L\(tmp_rhs)))';
    
elseif( 1 == 1 )
    b = uold(1,:)';
    up(2,inds) = c^2*(U\(L\(Kxx*b(inds))))';
else
    b = uold(1,:)';
    up(2,inds) = c^2*(M\(Kxx*b(inds)))';
    
end

if( closure ~= ghostBasis )
    nGhost = (p-1)/2;
    ileft  = 1+nGhost;
    iright = N-nGhost;
    if( iBCLeft == Neumann )
        for ig = 1:nGhost
            up(2,ileft-ig)  = up(2,ileft+ig);
        end
    else % iBCLeft == Dirichlet
        for ig = 1:nGhost
            up(2,ileft-ig)  = -up(2,ileft+ig);
        end
    end
    
    if( iBCRight == Neumann )
        for ig = 1:nGhost
            up(2,iright+ig) = up(2,iright-ig);
        end
    else % iBCRight == Dirichlet
        for ig = 1:nGhost
            up(2,iright+ig) = -up(2,iright-ig);
        end
    end
end

return
