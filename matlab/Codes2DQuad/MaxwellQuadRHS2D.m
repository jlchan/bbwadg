function [rhsHQ, rhsEQ] = MaxwellQuadRHS2D(HQ,EQ)

% function [rhsHQ, rhsEQ] = MaxwellRHS2D(HQ,EQ)
% Purpose  : Evaluate RHS flux in 2D Maxwell TM form 

Globals2D;

Hx = HQ(:,:,1);
Hy = HQ(:,:,2);
Ez = EQ(:,:);

% Define field differences at faces
dHx = zeros(Nfp*Nfaces,K); dHx(:) = Hx(vmapM)-Hx(vmapP);
dHy = zeros(Nfp*Nfaces,K); dHy(:) = Hy(vmapM)-Hy(vmapP); 
dEz = zeros(Nfp*Nfaces,K); dEz(:) = Ez(vmapM)-Ez(vmapP);

% Impose reflective boundary conditions (Ez+ = -Ez-)
dHx(mapB) = 0; dHy(mapB) = 0; dEz(mapB) = 2*Ez(vmapB);

% evaluate upwind fluxes
alpha = 1.0; 
ndotdH =  nx.*dHx+ny.*dHy;
fluxHx =  ny.*dEz + alpha*(ndotdH.*nx-dHx);
fluxHy = -nx.*dEz + alpha*(ndotdH.*ny-dHy);
fluxEz = -nx.*dHy + ny.*dHx - alpha*dEz;

% local derivatives of fields
[Ezx,Ezy] = GradQuad2D(Ez); 
[CuHx,CuHy,CuHz] = CurlQuad2D(Hx,Hy,[]);

% compute right hand sides of the PDE's
rhsHx = -Ezy  + LIFT*(Fscale.*fluxHx)/2.0;
rhsHy =  Ezx  + LIFT*(Fscale.*fluxHy)/2.0;
rhsEz =  CuHz + LIFT*(Fscale.*fluxEz)/2.0;

rhsHQ(:,:,1) = rhsHx;
rhsHQ(:,:,2) = rhsHy;
rhsEQ = rhsEz;

return;
