% N = order of approx
% VX, VY, VZ = global vertex pos
% v = vertex positions for the element
% fv = vertex face nodes for the element

function [rf sf tf wf fids] = surface_cubature(N,VX,VY,VZ,fv)

hybridgGlobalFlags % for useSEM flag

Nfaces = length(fv);
rf = []; sf = []; tf = []; wf = [];
fids = cell(Nfaces,1);
off = 0;
for f = 1:Nfaces
    vf = fv{f};
    vx = VX(vf); vy = VY(vf); vz = VZ(vf);
    if length(vf)==4
        % GQ weights/cubature
        [r1Dq w1Dr] = JacobiGQ(0,0,N);
        s1Dq = r1Dq; w1Ds = w1Dr;
        if useSEM
            % GLL weights/cubature
            r1Dq = JacobiGL(0,0,N); s1Dq = r1Dq;
            V1D = Vandermonde1D(N,r1Dq);  M1D = inv(V1D*V1D');
            w1Dr = sum(M1D,2); % get GLL weights by lumping
            w1Ds = w1Dr;
        end

        if nodalLIFT % replace r quadrature with exact quadrature            
            [r1Dq w1Dr] = JacobiGQ(0,0,N);
            if ~useSEM 
                [s1Dq w1Ds] = JacobiGQ(0,0,N); % exact quadrature in vertical direction
            end
        end
        [rqq sqq] = meshgrid(r1Dq,s1Dq); 
        rqq = rqq(:);  sqq = sqq(:);
        [wr ws] = meshgrid(w1Dr,w1Ds); 
        wr = wr(:);  ws = ws(:); wqf = wr.*ws; 
        [rqf sqf tqf] = map_quad(vx,vy,vz,rqq,sqq);        
    elseif length(vf)==3                
        [rqf sqf wqf] = Cubature2D(2*N+1); % for *exact* wedges
        [rqf sqf tqf] = map_triangle(vx,vy,vz,rqf,sqf);                
    end 
    wqf = wqf/sum(wqf); % take care of scaling in sJ
    
    fids{f} = off + (1:length(rqf));
    off = off + length(rqf);
    rf = [rf; rqf(:)]; sf = [sf; sqf(:)]; tf = [tf; tqf(:)]; 
    wf = [wf; wqf(:)];        
end