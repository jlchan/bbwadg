function z = dphi( xi,p )
  % derivate for basis function for p-th order reconstruction
  
  if( p == 7 & 1 == 0 )
    z = 0;
    
    if( xi > 4 )
      z = 0;
    elseif( xi > 3 )
      z = 0.86e2 / 0.2553e4 * xi - 0.301e3 / 0.2553e4 - 0.7e1 / 0.666e3 * (xi - 0.4e1) * (xi - 0.5e1) - 0.7e1 / 0.666e3 * (xi - 0.3e1) * (xi - 0.5e1) - 0.7e1 / 0.666e3 * (xi - 0.3e1) * (xi - 0.4e1);
    elseif( xi > 2 )
      z = 0.109e3 / 0.5106e4 * xi - 0.545e3 / 0.10212e5 - 0.5e1 / 0.1332e4 * (xi - 0.3e1) * (xi - 0.4e1) - 0.5e1 / 0.1332e4 * (xi - 0.2e1) * (xi - 0.4e1) - 0.5e1 / 0.1332e4 * (xi - 0.2e1) * (xi - 0.3e1);
    elseif( xi > 1 )
      z = 0.22e2 / 0.2553e4 * xi - 0.11e2 / 0.851e3 - (xi - 0.2e1) * (xi - 0.3e1) / 0.1332e4 - (xi - 0.1e1) * (xi - 0.3e1) / 0.1332e4 - (xi - 0.1e1) * (xi - 0.2e1) / 0.1332e4;
    elseif( xi > 0 )
      z = -0.1851e4 / 0.1702e4 + 0.149e3 / 0.851e3 * xi + 0.7e1 / 0.74e2 * (xi - 0.1e1) * (xi - 0.2e1) + 0.7e1 / 0.74e2 * xi * (xi - 0.2e1) + 0.7e1 / 0.74e2 * xi * (xi - 0.1e1);
    elseif( xi > -1 )
      z = 0.684e3 / 0.851e3 - 0.334e3 / 0.851e3 * xi - 0.7e1 / 0.74e2 * xi * (xi - 0.1e1) - 0.7e1 / 0.74e2 * (xi + 0.1e1) * (xi - 0.1e1) - 0.7e1 / 0.74e2 * (xi + 0.1e1) * xi;
    elseif( xi > -2 )
      z = 0.67e2 / 0.5106e4 * xi + 0.67e2 / 0.3404e4 + (xi + 0.1e1) * xi / 0.1332e4 + (xi + 0.2e1) * xi / 0.1332e4 + (xi + 0.2e1) * (xi + 0.1e1) / 0.1332e4;
    elseif( xi > -3 )
      z = 0.112e3 / 0.2553e4 * xi + 0.280e3 / 0.2553e4 + 0.5e1 / 0.1332e4 * (xi + 0.2e1) * (xi + 0.1e1) + 0.5e1 / 0.1332e4 * (xi + 0.3e1) * (xi + 0.1e1) + 0.5e1 / 0.1332e4 * (xi + 0.3e1) * (xi + 0.2e1);
    elseif( xi > -4 )
      z = 0.247e3 / 0.2553e4 * xi + 0.1729e4 / 0.5106e4 + 0.7e1 / 0.666e3 * (xi + 0.3e1) * (xi + 0.2e1) + 0.7e1 / 0.666e3 * (xi + 0.4e1) * (xi + 0.2e1) + 0.7e1 / 0.666e3 * (xi + 0.4e1) * (xi + 0.3e1);
    end
    return
    
  end
    
  if( p == 5 & 1 == 0 )
    z = 0;
    
    if( xi > 3 )
      z = 0;
    elseif( xi > 2 )
      z = 0.7e1 / 0.130e3 * xi - 0.7e1 / 0.52e2 - 0.5e1 / 0.156e3 * (xi - 0.3e1) * (xi - 0.4e1) - 0.5e1 / 0.156e3 * (xi - 0.2e1) * (xi - 0.4e1) - 0.5e1 / 0.156e3 * (xi - 0.2e1) * (xi - 0.3e1);
    elseif( xi > 1 )
      z = 0.2e1 / 0.65e2 * xi - 0.3e1 / 0.65e2 - (xi - 0.2e1) * (xi - 0.3e1) / 0.156e3 - (xi - 0.1e1) * (xi - 0.3e1) / 0.156e3 - (xi - 0.1e1) * (xi - 0.2e1) / 0.156e3;
    elseif( xi > 0 )
      z = -0.76e2 / 0.65e2 + 0.22e2 / 0.65e2 * xi + 0.7e1 / 0.39e2 * (xi - 0.1e1) * (xi - 0.2e1) + 0.7e1 / 0.39e2 * xi * (xi - 0.2e1) + 0.7e1 / 0.39e2 * xi * (xi - 0.1e1);
    elseif( xi > -1 )
      z = 0.41e2 / 0.65e2 - 0.48e2 / 0.65e2 * xi - 0.7e1 / 0.39e2 * xi * (xi - 0.1e1) - 0.7e1 / 0.39e2 * (xi + 0.1e1) * (xi - 0.1e1) - 0.7e1 / 0.39e2 * (xi + 0.1e1) * xi;
    elseif( xi > -2 )
      z = 0.9e1 / 0.130e3 * xi + 0.27e2 / 0.260e3 + (xi + 0.1e1) * xi / 0.156e3 + (xi + 0.2e1) * xi / 0.156e3 + (xi + 0.2e1) * (xi + 0.1e1) / 0.156e3;
    elseif( xi > -3 )
      z = 0.16e2 / 0.65e2 * xi + 0.8e1 / 0.13e2 + 0.5e1 / 0.156e3 * (xi + 0.2e1) * (xi + 0.1e1) + 0.5e1 / 0.156e3 * (xi + 0.3e1) * (xi + 0.1e1) + 0.5e1 / 0.156e3 * (xi + 0.3e1) * (xi + 0.2e1);
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
      num = 0;
      for k = 1:length(lp)
        den = den*lp(k);

        tmp = 1;
        for i = 1:length(lp)
          if( i ~= k )
            tmp = tmp*(xi+lp(i));
	  end
        end
        num = num+tmp;
      end
      z = num/den;
    end

  else
    z = 0;
    disp( 'error' );
    pause;
  end
