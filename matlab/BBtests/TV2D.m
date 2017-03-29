function [TV TVx TVy] = TV2D(u)

Globals2D

off = 0;
TVx = 0;
TVy = 0;
for i = 0:N
    Ni = N-i;
    Npi = Ni+1;
    idsx = (1:Npi) + off;
    TVx = TVx + TV1D(N-i,u(idsx,:));

    offj = 0;
    idsy = [];
    for j = 0:N-i
        idsy(j+1) = i + offj + 1;
        offj = offj + (N-j+1);
    end    
    TVy = TVy + TV1D(N-i,u(idsy,:));
    off = off + Npi;
end
TV = TVx+TVy;


function TV = TV1D(N,u)

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end