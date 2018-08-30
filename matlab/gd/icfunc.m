function z = icfunc( xi,icOption )

  if( icOption == 1 )
    alpha = 20;
    z = exp(-alpha*xi^2);
  elseif( icOption == 2 )
    if( xi < 0 & xi > -0.25 )
      z = xi+0.25;
    elseif( xi >= 0 & xi < 0.25 )
      z = -xi+0.25;
    else
      z = 0;
    end
  elseif( icOption == 3 )
    if( abs(xi) < 0.25 )
      z = 1;
    else
      z = 0;
    end
  else
    z = sin(pi*xi);
  end

  return
  end
