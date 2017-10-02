% builds maps between xf,yf,zf
function [mapM,mapP,mapB] = BuildHybridMaps()

hybridgGlobalFlags
hybridgGlobals3D

nodeids = reshape(1:K*NfcMax, NfcMax, K);
mapM   = nodeids; % initialize to node ids
mapP   = zeros(NfcMax, K);

for e=1:K
    
    if Nv(e)==4 % tet
        fids = fidsT;
    elseif Nv(e) == 5 % pyr
        fids = fidsP;
    elseif Nv(e) == 6 % wedge
        fids = fidsW;
    elseif Nv(e) == 8 % hex
        fids = fidsH;
    end
    
    Nfaces = length(fids);
    for f = 1:Nfaces        
        enbr = EToE(e,f);
        if (enbr==0) % boundary node
            mapP(nodeids(fids{f},e)) = mapM(nodeids(fids{f},e));
        else
            fnbr = EToF(e,f);
            
            if Nv(enbr)==4 % tet
                fnbrids = fidsT;
            elseif Nv(enbr) == 5 % pyr
                fnbrids = fidsP;
            elseif Nv(enbr) == 6 % wedge
                fnbrids = fidsW;
            elseif Nv(enbr) == 8 % hex
                fnbrids = fidsH;
            end
            
            idM = nodeids(fids{f},e); idM = idM(:);
            idP = nodeids(fnbrids{fnbr},enbr); idP = idP(:);
            
%             clf
%             plot3(xf(idM),yf(idM),zf(idM),'.');
%             hold on
%             plot3(xf(idP),yf(idP),zf(idP),'ro');
%             title(sprintf('elem %d, face %d',e,f))
%             pause
            
            tmp = ones(1,length(fids{f}));
            xM = xf(idM)*tmp; yM = yf(idM)*tmp; zM = zf(idM)*tmp;
            xP = xf(idP)*tmp; yP = yf(idP)*tmp; zP = zf(idP)*tmp;
            
            % Compute distance matrix            
            D = (xM -xP').^2 + (yM-yP').^2 + (zM-zP').^2;
            
            [idM, idP] = find(abs(D)<NODETOL);
            mapP(fnbrids{fnbr}(idP),enbr) = mapM(fids{f}(idM),e);
        end
        
    end
end
% mapP = mapP(:); mapM = mapM(:);

mapM = reshape(mapM,NfcMax,K);
mapP = reshape(mapP,NfcMax,K);

% Create list of boundary nodes
mapB = mapM(mapP==mapM);