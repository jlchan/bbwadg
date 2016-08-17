function rate = compute_rate(err,ids)
if nargin==1
    ids = 1:length(err);
end
h = .5.^(1:length(err)); 
h = h(:); 
err = err(:);
fit = [log(h(ids)) ones(size(h(ids)))]\log(err(ids));
%title(sprintf('Order %d: rate = %f\n',N,fit(1)))
rate = fit(1);
