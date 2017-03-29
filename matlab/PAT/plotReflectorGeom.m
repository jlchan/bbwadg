load PAT/TripleHair_Param
load PAT/TripleHair_Data
load PAT/TripleHair_Recon

xn = xn + x0_shift;
x0_img = x0_img + x0_shift;
% color_line3(x0_img(:),-y0_img(:),Reconn_p0,'b.','DisplayName','reconstruc
% tion zone')
vv = Reconn_p0;
color_line3(x0_img(:),y0_img(:),vv,vv,'.')
hold on
plot(xn,yt*ones(size(xn)),'r.','DisplayName','receivers')
text(mean(x0_img(:))+5,mean(y0_img(:)),sprintf('reconstruction \n zone'))
text(mean(xn),yt+2,'receivers')

plot(ym_left,0,'k*')
text(ym_left+2,-2,sprintf('left mirror \n "intercept"'))
plot(ym_right,0,'k*')
text(ym_right+2,-2,sprintf('right mirror \n "intercept"'))

yintercept=0;

% left mirror line y = ax+b, 
a = -1/tan(pi/2-alpha);
b = yintercept-a*ym_left; x = linspace(-75,25,1000);
plot(x,a*x+b,'k-')

% right mirror line y = ax+b, 
yintercept=ym_left;
a = 1/tan(pi/2-beta);
b = yintercept-a*ym_right; x = linspace(-25,75,1000);
plot(x,a*x+b,'k-')

a = 75;
axis([-a a -a/2 a/2])
plot(linspace(-a,a,25),yt*ones(1,25),'r--')

plot(14.9071,0,'bs')
plot(35.3,yt,'rs')
plot(-10,yt,'ks')

plot(14.9071,-13.2433,'bs')
plot(55.6929,13.2433,'rs')
plot(-34.9070,13.2433,'ks')

