function [r s t w] = hex_cubature(N)

hybridgGlobalFlags

[r1D w1D] = JacobiGQ(0,0,N);
if useSEM
    r1D = JacobiGL(0,0,N);
    V1D = Vandermonde1D(N,r1D);  M1D = inv(V1D*V1D');
    w1D = sum(M1D,2); % get GLL weights by lumping
end
[r s t] = meshgrid(r1D);
[wr ws wt] = meshgrid(w1D);
r = r(:); s = s(:); t = t(:);
wr = wr(:); ws = ws(:); wt = wt(:);
w = wr.*ws.*wt;