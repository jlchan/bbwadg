function test_bern_multiply

N = 3;
u = (1:N+1)';
%v = (0:N)';
v = zeros(N+1,1);
v(1) = 1;

cN = bern_coeffs(N);
uc = u.*cN; vc = v.*cN;

[r w] = JacobiGQ(0,0,2*N);
V = bern_basis_1D(N,r);
% Vs = bern_scaled(N,r);

V2N = bern_basis_1D(2*N,r);
uv = V2N\((V*u).*(V*v));

norm(V2N*uv-(V*u).*(V*v))

c2N = bern_coeffs(2*N);
% uvc = zeros(2*N+1,1);
% for j = 1:N+1    
%     uvc(j:j+N) = uvc(j:j+N) + v(j)*u;
% end
% uvc = uvc./c2N;

uvc= conv(uc,vc)./c2N;


M = V'*diag(w)*V;
I = M\(V'*diag(w)*V2N);
uvt = I*uvc;
% I = I./min(abs(I(:)));
% I(1:N+1,1:N+1)
rp = linspace(-1,1,250);
Vp = bern_basis_1D(N,rp);
V2Np = bern_basis_1D(2*N,rp);

plot(rp,V2Np*uv,'b-');hold on
plot(rp,Vp*uvt,'r-')
%uvt = ifft(fft(uc,2*N+1).*fft(vc,2*N+1))./c2N;
plot(rp,Vp*uvt,'ko-')

keyboard

function V = bern_scaled(N,r)

r = (1+r)/2; % convert to unit
for i = 0:N
    V(:,i+1) = (r.^i).*(1-r).^(N-i);  
end

function c = bern_coeffs(N)

for i = 0:N
    c(i+1,1) = nchoosek(N,i);
end
