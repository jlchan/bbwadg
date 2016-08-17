function plotConvergence(h,err,pos,expectedOrder)

fit = [log(h(:)) ones(size(h(:)))]\log(err(:));
order = fit(1);
if nargin<3
    pos = 'top';
end
if nargin==4
    order = expectedOrder;
end

if strcmp(pos,'top')
%     % Draw the solution order
%     text(h(1)+5.e-2,err(1),sprintf('N = %d',N),'FontSize',14)

    x1 = h(end); x2 = h(end-1)-0.5*(h(end-1)-h(end));
    y1 = err(end) + err(end)*0.2;
    b = log(y1) - order*log(x1);
    y2 = exp(order*log(x2) + b);
    x3 = x1; y3 = y2;
    px = [x1,x2,x3,x1];py = [y1,y2,y3,y1];
    
    % plot the order of convergence
    xm = x1; xm = xm - 0.2*xm;
    ym = y1 + 0.3*(y3-y1);
    text(xm*.85,ym,num2str(order,3),'FontSize',14)
    
    xm = x1 + 0.3*(x2-x1);
    ym = y2 + 0.5*y2;
%     text(xm,ym,'1','FontSize',14)
    loglog(px,py,'b','linewidth',1.0)

else
    x1 = h(end); x2 = h(end-1)-0.5*(h(end-1)-h(end));
    y1 = .8*err(end);
    b = log(y1) - order*log(x1);
    y2 = exp(order*log(x2) + b);
    x3 = x2; y3 = y1;
    px = [x1,x2,x3,x1];py = [y1,y2,y3,y1];
    
    % plot the order of convergence
    xm = x2; xm = .8*xm;
    ym = .65*y1;
    text(xm*.85,ym,num2str(order,3),'FontSize',14)
    
    xm = x2 + 0.2*(x2-x1);
    ym = y2 + 0.75*(y3-y2);
%     text(xm,ym,'1','FontSize',14)
    loglog(px,py,'r','linewidth',1.0)

end
%text(px,py,num2str((1:4)'))
% keyboard



return

