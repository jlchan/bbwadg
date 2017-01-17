function [rf sf tf wf fids] = wedge_surface_cubature(N)

% SEM nodes/weights - lumped cubature
% [t1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,t1D); w1D = sum(inv(V*V'),2);
[t1D w1D] = JacobiGQ(0,0,N);  % DENSE LIFT MATRICES

% face quadrature nodes
off = 0;

% tri faces (face 1/5)
[rff sff wff] = Cubature2D(2*N+1); tff = -ones(size(rff)); 
fids = 1:length(rff);
rf(off+fids) = rff; 
sf(off+fids) = sff; 
tf(off+fids) = tff; 
wf(off+fids) = wff;
ffids{1} = off + fids;
off = off + length(rff);

% quad faces: TP of GLL and whatever in r/s direction 
[s1D w1Ds] = JacobiGQ(0,0,N); 
% [s1D] = JacobiGL(0,0,N); w1Ds = w1D; % GLL quadrature 
% [s1D] = JacobiGQ(0,0,N); w1Ds = w1D; % GQ nodes/GLL quadrature - why the hell does this work so well?

% face s = -1
[rff tff] = meshgrid(s1D,t1D); rff = rff(:); tff = tff(:); 
sff = -ones(size(rff)); [wrf wtf] = meshgrid(w1Ds,w1D); wff = wrf(:).*wtf(:);
fids = 1:length(rff);
rf(off+fids) = rff; 
sf(off+fids) = sff; 
tf(off+fids) = tff; 
wf(off+fids) = wff;
ffids{2} = off + fids;
off = off + length(rff);

% face r+s = 0
[rff tff] = meshgrid(s1D,t1D); rff = rff(:); tff = tff(:); 
sff = -rff; [wrf wtf] = meshgrid(w1Ds,w1D); wff = wrf(:).*wtf(:); % don't scale - sJ takes care
fids = 1:length(rff); 
rf(off+fids) = rff; 
sf(off+fids) = sff;
tf(off+fids) = tff; 
wf(off+fids) = wff;
ffids{3} = off + fids;
off = off + length(rff);

% face r = -1
[sff tff] = meshgrid(s1D,t1D); sff = sff(:); tff = tff(:);  
rff = -ones(size(sff)); [wrf wtf] = meshgrid(w1Ds,w1D); wff = wrf(:).*wtf(:);
fids = 1:length(rff);
rf(off+fids) = rff; 
sf(off+fids) = sff; 
tf(off+fids) = tff; 
wf(off+fids) = wff;
ffids{4} = off + fids;
off = off + length(rff);

% face 5: t = 1
[rff sff wff] = Cubature2D(2*N+1); tff = ones(size(rff));
fids = 1:length(rff);
rf(off+fids) = rff; 
sf(off+fids) = sff; 
tf(off+fids) = tff; 
wf(off+fids) = wff;
ffids{5} = off + fids;
off = off + length(rff);

fids = ffids; 
