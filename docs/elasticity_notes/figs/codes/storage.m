clear
K = 500000;
for N = 1:6
   Nptri = (N+1)*(N+2)/2;
   storageLIFT(N) = Nptri^2 + Nptri*(N+1)*3;
   storageU(N) = 4*Nptri*(N+1); % 4 fields and N+1 triangles   
%    storageU(N) = 4*(N+1)^3; % 4 fields and TP
end
plot(storageU*K/1e9,'bo-','linewidth',2); hold on; plot(storageLIFT*K/1e9,'rs-','linewidth',2)
% plot(storageU./storageLIFT,'bo-','linewidth',2);