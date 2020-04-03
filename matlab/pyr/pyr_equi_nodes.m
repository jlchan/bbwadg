function [r s t] = pyr_equi_nodes(N)

r1D = linspace(-1,1,N+1)';
r = []; s = []; t = [];
for level = 0:N
    rlevel = r1D(1:N+1-level);
    [rr ss] = meshgrid(rlevel);
    r = [r; rr(:)];
    s = [s; ss(:)];    
    t = [t; r1D(level+1)*ones(size(rr(:)))];
end

