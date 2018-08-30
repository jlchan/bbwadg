function c = generalDifferenceWeights( z,x,m ) 

c1 = 1;
c4 = x(1)-z;

n = length(x)-1;
c = zeros(n+1,m+1);

c(1,1) = 1;
for i = 1:n 
  mn = min(i,m);
  c2 = 1;
  c5 = c4;

  c4 = x(i+1)-z;
  for j =0:i-1
    c3 = x(i+1)-x(j+1);
    c2 = c2*c3;
    if (j == (i-1) ) 
      for k = mn:-1:1
        c(i+1,k+1) = c1*(k*c(i-1+1,k-1+1)-c5*c(i-1+1,k+1))/c2;
      end
      c(i+1,0+1) = -c1*c5*c(i-1+1,1)/c2;
    end

    for k = mn:-1:1
      c(j+1,k+1) = (c4*c(j+1,k+1)-k*c(j+1,k-1+1))/c3;
    end
    c(j+1,1) = c4*c(j+1,1)/c3;
  end
  c1 = c2; 
end

return
end