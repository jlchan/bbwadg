% 2D

clear
clc
N = 3;
Np = (N+1)*(N+2)/2;
id = [];
jjsk = 0;
for jj = 0:N
    for ii = 0:N-jj
        jsk = 0;
        for j = 0:N
            for i = 0:N-j
                id = [id (i + ii + jsk + jjsk + 1)];
            end
            jsk = jsk + (2*N-j-jj+1);
        end
    end
    jjsk = jjsk + (2*N-jj+1);
end

id

[r s] = EquiNodes2D(2*N); [r s] = xytors(r,s);
plot(r,s,'o')
text(r+.1,s,num2str((1:length(r))'))

%% 3d

clear
clc
N = 2;
Np = (N+1)*(N+2)*(N+3)/6;
id = [];
kksk = 0;
for kk = 0:N
    jjsk = 0;
    for jj = 0:N-kk
        for ii = 0:N-jj-kk
            
            ksk = 0;
            for k = 0:N
                jsk = 0;
                for j = 0:N-k
                    for i = 0:N-j-k
                        idtmp = (i + ii + jsk + jjsk + ksk + kksk + 1);
                        id = [id  idtmp];
                        [idtmp i jsk ksk jjsk kksk]
                    end
                    NN = 2*N- j-jj - k-kk;
                    jsk = jsk + (NN+1);
                end
                NN = 2*N - k-kk;
                ksk = ksk + (NN+1)*(NN+2)/2 - jj;
            end
            
        end
        NN = 2*N-jj-kk;
        jjsk = jjsk + (NN+1);
    end
    NN = 2*N - kk;
    kksk = kksk + (NN+1)*(NN+2)/2;
end

id

[r s t] = EquiNodes3D(2*N);
for i = 1:length(id)
    clf
    plot3(r,s,t,'o')
    hold on;plot3(r(id(i)),s(id(i)),t(id(i)),'x')
    text(r+.05,s,t,num2str((1:length(r))'))
    view(5,10);%return
    pause
end