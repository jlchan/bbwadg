function [u Klim alpha] = limit1D(u)

[u Klim alpha] = p1limit(u);


function [ulimit ids alpha] = p1limit(u)

Globals1D
global VB

% Compute cell averages
uh = invV*u;
uh(2:Np,:)=0;
uavg = V*uh;
v = uavg(1,:);

% Apply slope limiter as needed.
ulimit = u;
eps0=1.0e-8;
% find end values of each element
ue1 = u(1,:);
ue2 = u(end,:);

% find cell averages
vk = v;
%vkm1 = [v(1),v(1:K-1)]; vkp1 = [v(2:K),v(K)];
vkm1 = [v(K),v(1:K-1)]; vkp1 = [v(2:K),v(1)]; % periodic

% Apply reconstruction to find elements in need of limiting
ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

if 0
    % convert to bernstein, find K with large TV
    uB = VB\u;
    TV = 0;
    for i = 1:N
        TV = TV + abs(uB(i,:) - uB(i+1,:));
    end
    TV = TV./(N*max(abs(uB(:))));
    ids = find(TV > .5*max(TV));
end

% Check to see if any elements require limiting
if(~isempty(ids))
    % create piecewise linear solution for limiting on specified elements
    uhl = invV*u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
    % apply slope limiter to selected elements
    ulimit(:,ids) = SlopeLimitLin(ul,x(:,ids),vkm1(ids),vk(ids),vkp1(ids));
end
alpha = 0*ids;

return;

function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1)
% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right
Globals1D;

% Compute various geometric measures
% ulimit = ul;
h = xl(Np,:)-xl(1,:);
x0 = ones(Np,1)*(xl(1,:) + h/2);
hN = ones(Np,1)*h;

% Limit function
ux = (2./hN).*(Dr*ul);
ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*minmod([ux(1,:); (vp1-v0)./h; (v0-vm1)./h]));

return


function mfunc = minmod(v)
% function mfunc = minmodB(v,M,h)
% Purpose: Implement the TVB modified midmod function. v is a vector

Globals1D
M = 50; h = VX(2)-VX(1);
mfunc = v(1,:);
ids = find(abs(mfunc) > M*h.^2);
if(size(ids,2)>0)
    mfunc(ids) = minmod_orig(v(:,ids));
end
return


function mfunc = minmod_orig(v)
% function mfunc = minmod(v)
% Purpose: Implement the midmod function v is a vector
m = size(v,1);
mfunc = zeros(1,size(v,2));
s = sum(sign(v),1)/m;
ids = find(abs(s)==1);
if(~isempty(ids))
    mfunc(ids) = s(ids).*min(abs(v(:,ids)),[],1);
end
return;