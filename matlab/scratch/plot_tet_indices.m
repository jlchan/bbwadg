clear
N = 2;
[r s t] = EquiNodes3D(N);  %[r s t] = xyztorst(r,s,t);

% find all the nodes that lie on each edge
NODETOL = 1e-6;
fmask1   = find( abs(1+t) < NODETOL)'; fmask2   = find( abs(1+s) < NODETOL)';
fmask3   = find( abs(1+r+s+t) < NODETOL)'; fmask4   = find( abs(1+r) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

sk = 1;
for l = 0:N
    for k = 0:N-l
        for j = 0:N-l-k
            i = N-j-k-l;
            idi(sk) = i;
            idj(sk) = j;
            idk(sk) = k;
            idl(sk) = l;
            sk = sk + 1;
        end
    end
end

plot3(r,s,t,'o')
hold on
for ii = 1:length(r)
    l = idl(ii); k = idk(ii); j = idj(ii);
%     id1 = (N+1)*(N+2)/2*l - (2*N+3)*(l-1)*(l)/4 + (l-1)*(l)*(2*l-1)/12;
%     id2 = (N+1-l)*k - (k-1)*k/2;
%     id3 = j;
%     vol_id = id1+id2+id3+1;
    map3D = @(N,j,k,l) (N+1)*(N+2)/2*l - (2*N+3)*(l-1)*(l)/4 + (l-1)*(l)*(2*l-1)/12 + (N+1-l)*k - (k-1)*k/2 + j + 1;
    text(r(ii)+.01,s(ii),t(ii),sprintf('(%d,%d,%d,%d,%d,%d)',l,k,j,l,map3D(N,j,k,l),ii))
end
%%
f = 1;
rf = r(Fmask(:,f));
sf = s(Fmask(:,f));
tf = t(Fmask(:,f));
%plot(rf,sf,'r*')
plot3(rf,sf,tf,'r*') 
for i = 1:length(rf)
    %text(rf(i)+.01,sf(i),sprintf('(%d,%d,%d,%d)',idi(Fmask(i,f)),idj(Fmask(i,f)),idk(Fmask(i,f)),idl(Fmask(i,f))))
    text(rf(i)+.01,sf(i),tf(i),sprintf('(%d,%d,%d,%d)',idi(Fmask(i,f)),idj(Fmask(i,f)),idk(Fmask(i,f)),idl(Fmask(i,f))))
end

%%
[r s] = EquiNodes2D(N); [r s] = xytors(r,s);

sk = 1;
for k = 0:N
    for j = 0:N-k
        i = N-j-k;
        idi(sk) = i;
        idj(sk) = j;
        idk(sk) = k;
        sk = sk + 1;
    end
end

plot(r,s,'o')
hold on
for i = 1:length(r)
    text(r(i)+.01,s(i),sprintf('(%d,%d,%d)',idi(i),idj(i),idk(i)))
end