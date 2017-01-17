% GQ cubature on quad
function [r,s,w] = QCubature2D(N)

[r, w] = JacobiGQ(0, 0, N);
r = r*ones(1,N+1);
wr = w*ones(1,N+1);
s = r';
ws = wr';
r = r(:); s = s(:);w = wr(:).*ws(:);
