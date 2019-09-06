% check lagrange formula for GD

p = 3;
q = (p+1)/2;

cj = zeros(2*q,1);
sk = 1;
for j = -q:-1
    cj(sk) = (-1).^q*(-1)^j*nchoosek(p,q+j)/factorial(p);
    cj(2*q-sk+1) = -cj(sk);
    sk = sk + 1;
end

sk = 1;
for j = -q:q-1
    
    x = linspace(j,j+1,100);
    prod = ones(size(x));
    for k = j-q+1:j+q
        if k~=0
            prod = prod.*(x-k);
        end        
    end
    plot(x,cj(sk)*prod,'linewidth',2)
    hold on
    sk = sk + 1;
end

cj