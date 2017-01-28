function [mapM, mapP] = BuildNodeMaps2D(xf,yf,EToE,EToF)

K = size(EToE,1);
Nfaces = 3;
NODETOL = 1e-5;

% number volume nodes consecutively
Nfp = size(xf,1) / Nfaces; 
mapM    = (1:K*Nfp*Nfaces)';     
mapP    = reshape(mapM, Nfp, Nfaces, K);

one = ones(1, Nfp);
for k1=1:K
    for f1=1:Nfaces
        % find neighbor
        k2 = EToE(k1,f1); f2 = EToF(k1,f1);
        
        if (k1==k2)
            mapP(idM,f1,k1) = mapM(idM,f1,k1);
        else
            
            % find find volume node numbers of left and right nodes
            idM = (1:Nfp) + (f1-1)*Nfp;
            idP = (1:Nfp) + (f2-1)*Nfp;
            x1 = xf(idM,k1); y1 = yf(idM,k1);
            x2 = xf(idP,k2); y2 = yf(idP,k2);
            x1 = x1*one;  y1 = y1*one;
            x2 = x2*one;  y2 = y2*one;
            
            % Compute distance matrix
            D = (x1 -x2').^2 + (y1-y2').^2;
            [idM, idP] = find(sqrt(abs(D)) < NODETOL);
            mapP(idM,f1,k1) = idP + (f2-1)*Nfp+(k2-1)*Nfaces*Nfp;
        end
    end
end

mapP = mapP(:);
return
