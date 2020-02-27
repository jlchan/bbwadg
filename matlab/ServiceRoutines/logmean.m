% computes (aL-aR) / (log(aL)-log(aR))

function val = logmean(aL,aR,logL,logR)

if 0
    xi = aL./aR;
    f = (xi-1)./(xi+1);
    u = f.^2;
    ids = abs(u) < 1e-4; % arbitrary
    ui = u(ids);
    uu = ui.*ui;
    uuu = ui.*uu;
    uuuu = ui.*uuu;
    F = zeros(size(aL));
    F(~ids) = log(xi(~ids))./2./f(~ids);
    F(ids) = 1 + ui*.333333333333 + uu*.2 + uuu*0.142857142857143 + uuuu*.11111111111111;
    val = (aL+aR)./(2*F);
    
else
    
    % andrew winters approach
    
    da = aR-aL;
    aavg = .5*(aR+aL);
    f = da./aavg;
    v = f.^2;
    ids = abs(f) < 1e-4;
    vv = v(ids);
    val = zeros(size(aL));
    if nargin==4        
        val(~ids) = (aL(~ids)-aR(~ids))./(logL(~ids)-logR(~ids));
    else
        val(~ids) = (aL(~ids)-aR(~ids))./(log(aL(~ids))-log(aR(~ids)));
    end
    val(ids) = aavg(ids).*(1 + vv.*(-.2-vv.*(.0512 - vv*0.026038857142857)));
%     v2 = v(ids).*v(ids);
%     val(ids) = aavg(ids).*(1 - v(ids).*.2-v2.*.0512 + v2.*v(ids)*0.026038857142857);
end

% gamma = 1.4;
% c1 = (gamma-2)/3;
% c2 = (gamma+1)*(gamma-2)*(gamma-3)/45;
% c3 = (gamma+1)*(gamma-2)*(gamma-3)*(2*gamma*(gamma-2)-9)/945;
    




