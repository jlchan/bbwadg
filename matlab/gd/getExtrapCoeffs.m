function c = getExtrapCoeffs( order )
  if( order == 0 )
    c = [1];
  elseif( order == 1 )
    c = [1,1];
  else
    cold = [1,1];
    for i = 2:order
      cnew = [1];
      for k = 2:i
        cnew = [cnew,cold(k-1)+cold(k)];
      end
      cnew = [cnew,1];
      cold = cnew;
    end
    c = cold;
  end
  
  for k = 2:length(c)
    c(k) = (-1)^k*c(k);
  end
  
  return
  end