% scratchpad for checking ideas on BB basis

N = 3;

% [rq sq tq w] = tet_cubature(2*N);  
% Vq = bern_basis_tet(N,rq,sq,tq);

[r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
[rp sp tp] = EquiNodes3D(50);

[Vp] = bern_basis_tet(N,rp,sp,tp);
[V Vr Vs Vt V1 V2 V3 V4 id] = bern_basis_tet(N,r,s,t);

Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
Dt = V\Vt; Dt(abs(Dt)<1e-8) = 0;
D1 = V\V1; D1(abs(D1)<1e-8) = 0;
D2 = V\V2; D2(abs(D2)<1e-8) = 0;
D3 = V\V3; D3(abs(D3)<1e-8) = 0;
D4 = V\V4; D4(abs(D4)<1e-8) = 0;

sk = 1;
for i = 0:N
    for j = 0:N-i
        for k = 0:N-i-j
            l = N-i-j-k;
            idi(sk) = i;
            idj(sk) = j;
            idk(sk) = k;
            idl(sk) = l;
            sk = sk + 1;
        end
    end
end
idi = idi(id);
idj = idj(id);
idk = idk(id);
idl = idl(id);
ids = [idi;idj;idk;idl]

return
ids = [];
sk = 1;
for l = 0:N
    for k = 0:N-l
        for j = 0:N-l-k
            i = N-j-k-l;
            ids = [ids [i;j;k;l;sk]];
            sk = sk + 1;
        end
    end
end
ids

% reverse lookup
diff = 0;
for col = 1:length(id)
    i = ids(1,col);
    j = ids(2,col);
    k = ids(3,col);
    l = ids(4,col);
    
    off1 = 0;
    for ll = 0:l-1
        off1 = off1 + (N-ll+1)*(N-ll+2)/2;
    end
    off2 = 0;
    for kk = 0:k-1
        off2 = off2 + (N-l-kk+1);
    end
    
    err = abs((col-1)-(off1+off2+j));
    diff = diff + err;
    if err>0
        keyboard
    end
end
diff


return

color_line3(rp,sp,tp,Vp*Dr*(V\r),'.') % check drdr = 1


for i = 1:size(V,2);
    clf
    color_line3(rp,sp,tp,Vp(:,i),'.');hold on
    plot3(re(i),se(i),te(i),'o','markersize',20);
    view(3)
    pause
end
