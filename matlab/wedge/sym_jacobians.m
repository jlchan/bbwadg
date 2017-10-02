syms a b c
r = (1+a)*(1-c)/2 - 1;
s = b; 
t = c;
dadt = (2*(r + 1))/(t - 1)^2

%%
syms r s t

a = 2*(r+1)/(1-t) - 1;
b = s;
c = t;

J1 = [diff(a,r),diff(a,s),diff(a,t);
 diff(b,r),diff(b,s),diff(b,t);   
 diff(c,r),diff(c,s),diff(c,t)]

simplify(det(J1))

syms a b c
r = (a+1)*(1-c)/2 - 1;
s = b;
t = c;

J2 = [diff(r,a),diff(r,b),diff(r,c);
 diff(s,a),diff(s,b),diff(s,c);   
 diff(t,a),diff(t,b),diff(t,c)]

simplify(det(J2))

%% mapping for wedges

syms r s t

L1 = (1-s)/2.*(-(r+t))/2;
L2 = (1-s)/2.*(1+r)/2;
L3 = (1+s)/2.*(1+r)/2;
L4 = (1+s)/2.*(-(r+t))/2;
L5 = (1-s)/2.*(1+t)/2;
L6 = (1+s)/2.*(1+t)/2;

syms A B C D E F
r1 = [-A; 1; 1; -1; -1; -1]; s1 = [-B; -1; 1; 1; -1; 1]; t1 = [-C; -1; -1; -1; 1; 1];
 
L = [L1 L2 L3 L4 L5 L6];
xr = diff(L*r1,r);
xs = diff(L*r1,s);
xt = diff(L*r1,t);

yr = diff(L*s1,r);
ys = diff(L*s1,s);
yt = diff(L*s1,t);

zr = diff(L*t1,r);
zs = diff(L*t1,s);
zt = diff(L*t1,t);

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
pretty(simplify(expand(J*4)))

%%
syms a b c

simplify(subs(subs(subs(J,r,(a+1)*(1-c)/2 - 1),s,b),t,c))



