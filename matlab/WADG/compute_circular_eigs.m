clear lamEx
load WADG_eigs
x0 = sort(imag(lamWADG));
x0 = x0(x0 > -1e-14);
x0 = uniquetol(x0);
% x0 = x0(1:50);
% x0 = linspace(0,10,100);
% x = linspace(0,10,200);
% plot(x,besselj(n,x))
% hold on
for n = 0:5 % frequency    
    for i = 1:length(x0)       
        lamEx(i,1+n) = fzero(@(x)besselj(n,x),x0(i),optimset('TolX',1e-9));
    end
end
tol = 1e-8;
lamEx = lamEx(lamEx(:) > -tol);
lamEx = uniquetol(lamEx(:),tol);
% plot(0*lam(:),lam(:),'o')
plot(lamEx,'o')
hold on
plot(sort(x0),'x')
% lamDGh = lamDGi(lamDGi(:)>1e-4);