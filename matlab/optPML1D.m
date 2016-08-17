clear
clc
smax = linspace(.1,200,10);
for delta = .125;
    for p0sigma = 0:1
        for p0approx = 0%:p0sigma
            for i = 1:length(smax)
                l2err(i) = min(Wave1D(smax(i),delta,p0sigma,p0approx),1e1);
                i
            end
            l2err
            semilogy(smax,l2err,'o-','DisplayName',sprintf('delta = %f, p0sigma = %d, p0approx = %d',delta,p0sigma,p0approx));hold on            
        end
    end
end
legend show