shockVortexTp3_EC;

Nplot = round(sqrt(size(x,1)));

hold on

for e = 1:size(x,2)
    xe = reshape(x(:,e),Nplot,Nplot);
    ye = reshape(y(:,e),Nplot,Nplot);
    rhoe = reshape(rho(:,e),Nplot,Nplot);
    surface(xe,ye,rhoe)
    if (mod(e,1000)==0)
        fprintf('e = %d\n',e);
    end        
end

shading interp
axis off
axis equal
axis tight
colormap(gray)
colorbar
set(gca,'fontsize',15)
caxis([.8640 1.6033]) % t = .3
% caxis([.8680 1.3571]) % t = .7
% print(gcf,'-dpng','shockVortexTp3_EC.png')
% print(gcf,'-dpng','shockVortexTp7_EC.png')

% t = .3
% 0.8643    1.6033 % LF
% 0.8640    1.6012 % matrix diss 

% t = .7
% 0.9420    1.3571 % LF caxis
% 0.8680    1.3128 % matrix diss caxis
