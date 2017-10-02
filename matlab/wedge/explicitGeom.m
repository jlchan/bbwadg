[a b c] = meshgrid(linspace(-1,1,30));
a = a(:); b = b(:); c = c(:);
[r s t] = wedge_abctorst(a,b,c);

L1 = (1-s)/2.*(-(r+t))/2;
L2 = (1-s)/2.*(1+r)/2;
L3 = (1+s)/2.*(1+r)/2;
L4 = (1+s)/2.*(-(r+t))/2;
L5 = (1-s)/2.*(1+t)/2;
L6 = (1+s)/2.*(1+t)/2;

dL1r = s/4 - 1/4;
dL1s = r/4 + t/4;
dL1t = s/4 - 1/4;

dL2r = 1/4 - s/4;
dL2s = - r/4 - 1/4;
dL2t = 0*r;

dL3r = s/4 + 1/4;
dL3s = r/4 + 1/4;
dL3t = 0*r;

dL4r = - s/4 - 1/4;
dL4s = - r/4 - t/4;
dL4t = - s/4 - 1/4;

dL5r = 0*r;
dL5s = - t/4 - 1/4;
dL5t = 1/4 - s/4;

dL6r = 0*r;
dL6s = t/4 + 1/4;
dL6t = s/4 + 1/4;

x = [-1 1 1 -1 -1 -1]'; y = [-2 -1 1 1 -1 1]'; z = [-2 -1 -1 -1 1 1]';
V = [L1 L2 L3 L4 L5 L6];
Dr = [dL1r dL2r dL3r dL4r dL5r dL6r];
Ds = [dL1s dL2s dL3s dL4s dL5s dL6s];
Dt = [dL1t dL2t dL3t dL4t dL5t dL6t];


xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
color_line3(r,s,t,J,'.');colorbar;view(3)