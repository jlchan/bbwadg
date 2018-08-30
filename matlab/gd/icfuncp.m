function z = icfuncp( xi,icOption )

  if( icOption == 1 )
    alpha = 20;
    z = -2.0*alpha*xi*exp(-alpha*xi^2);
  elseif( icOption == 2 )
    if( xi < 0 & xi > -0.25 )
      z = 1;
    elseif( xi >= 0 & xi < 0.25 )
      z = -1;
    else
      z = 0;
    end
  elseif( icOption == 3 )
     z = 0;
  else
    z = pi*cos(pi*xi);
  end

  return
  end
