function z = phi( xi,p )
% basis function evaluation for p-th order reconstruction

if( p == 7 & 1 == 0 )
    z = 0;
    
    if( xi > 4 )
        z = 0;
    elseif( xi > 3 )
        z = 0.43e2 / 0.2553e4 * (xi - 3) * (xi - 4) - 0.7e1 / 0.666e3 * (xi - 3) * (xi - 4) * (xi - 5);
    elseif( xi > 2 )
        z = 0.109e3 / 0.10212e5 * (xi - 2) * (xi - 3) - 0.5e1 / 0.1332e4 * (xi - 2) * (xi - 3) * (xi - 4);
    elseif( xi > 1 )
        z = 0.11e2 / 0.2553e4 * (xi - 1) * (xi - 2) - ((xi - 1) * (xi - 2) * (xi - 3)) / 0.1332e4;
    elseif( xi > 0 )
        z = 0.1e1 - xi + 0.149e3 / 0.1702e4 * xi * (xi - 1) + 0.7e1 / 0.74e2 * xi * (xi - 1) * (xi - 2);
    elseif( xi > -1 )
        z = xi + 0.1e1 - 0.167e3 / 0.851e3 * (xi + 1) * xi - 0.7e1 / 0.74e2 * (xi + 1) * xi * (xi - 1);
    elseif( xi > -2 )
        z = 0.67e2 / 0.10212e5 * (xi + 2) * (xi + 1) + ((xi + 2) * (xi + 1) * xi) / 0.1332e4;
    elseif( xi > -3 )
        z = 0.56e2 / 0.2553e4 * (xi + 3) * (xi + 2) + 0.5e1 / 0.1332e4 * (xi + 3) * (xi + 2) * (xi + 1);
    elseif( xi > -4 )
        z = 0.247e3 / 0.5106e4 * (xi + 4) * (xi + 3) + 0.7e1 / 0.666e3 * (xi + 4) * (xi + 3) * (xi + 2);
    end
    return
    
end

if( p == 5 & 1 == 0 )
    z = 0;
    
    if( xi > 3 )
        z = 0;
    elseif( xi > 2 )
        z = .7e1 / 0.260e3 * (xi - 2) * (xi - 3) - 0.5e1 / 0.156e3 * (xi - 2) * (xi - 3) * (xi - 4);
    elseif( xi > 1 )
        z = ((xi - 1) * (xi - 2)) / 0.65e2 - ((xi - 1) * (xi - 2) * (xi - 3)) / 0.156e3;
    elseif( xi > 0 )
        z = 0.1e1 - xi + 0.11e2 / 0.65e2 * xi * (xi - 1) + 0.7e1 / 0.39e2 * xi * (xi - 1) * (xi - 2);
    elseif( xi > -1 )
        z = xi + 0.1e1 - 0.24e2 / 0.65e2 * (xi + 1) * xi - 0.7e1 / 0.39e2 * (xi + 1) * xi * (xi - 1);
    elseif( xi > -2 )
        z = 0.9e1 / 0.260e3 * (xi + 2) * (xi + 1) + ((xi + 2) * (xi + 1) * xi) / 0.156e3;
    elseif( xi > -3 )
        z = 0.8e1 / 0.65e2 * (xi + 3) * (xi + 2) + 0.5e1 / 0.156e3 * (xi + 3) * (xi + 2) * (xi + 1);
    end
    return
    
end

if( mod(p,2) == 1 )
    z = 0;
    a = -(p+1)/2;
    b = -a;
    if( xi >= a  & xi < b )
        ia = p-ceil(xi-a)+1;
        ib = -(p-ceil(b-xi)+1);
        lp = [ia:-1:1,-1:-1:ib];
        
        den = 1;
        num = 1;
        for k = 1:length(lp)
            num = num*(xi+lp(k));
            den = den*lp(k);
        end
        z = num/den;
    end
    
else
    disp( 'error' );
    pause;
end


