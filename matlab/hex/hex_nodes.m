function [r s t] = hex_nodes(N)
hybridgGlobalFlags
if useSEM
    r1D = JacobiGL(0,0,N);
%     [r s t] = meshgrid(JacobiGL(0,0,N));
else
    r1D = JacobiGQ(0,0,N);
%     [r s t] = meshgrid(JacobiGQ(0,0,N));
end

sk=1;
for i = 1:N+1
    for j = 1:N+1
        for k = 1:N+1
            r(sk) = r1D(k);
            s(sk) = r1D(j);
            t(sk) = r1D(i);
            sk = sk + 1;
        end
    end
end

r = r(:); s = s(:); t = t(:);