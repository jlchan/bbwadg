ref = 3;
surface_mesh = sprintf('Grid/tri%d.neu',ref);

%filename=['wedge_perturbxyz_' num2str(ref)];
% filename='ref_layer';
usurf2 = @(x,y) ones(size(x)); % assumes wavy layer is centered around z = 0.
usurf = @(x,y) -1*ones(size(x)); %.25*sin(2*x+2.5).*sin(3*y+3.5);
% nlevel = 1;

% filename='wavy_layer';
% usurf = @(x,y) .25*(sin(2*x+2.5).*sin(3*y+3.5) + .25*sin(4*x+1.5).*sin(5*y+5.5) + .125*sin(8*x).*sin(8*y));
% usurf2 = @(x,y) 1+.1*sin(x+2.5).*sin(3*y+3.5) + .05*cos(10*x).*sin(9*y);
nlevel = 2.^(ref-1);

height = 2;
[VXW VYW VZW EToVW EToEW EToFW KW num_interface_elems] = smoothed_wedge_layers(filename,surface_mesh,usurf,usurf2,nlevel);

if 1    
    h = (.5).^ref;
    
    % randomly perturb z-coords
    ids = abs(abs(VZW)-1) > 1e-8;
    VZW(ids) = VZW(ids) + h/2*randn(size(VZW(ids)));
    ids = abs(abs(VXW)-1) > 1e-8;
    VXW(ids) = VXW(ids) + h/2*randn(size(VXW(ids)));
    ids = abs(abs(VYW)-1) > 1e-8;
    VYW(ids) = VYW(ids) + h/2*randn(size(VYW(ids)));
    
%     t = linspace(-1,1,nlevel+1);
%     for lev = 1:nlevel+1
%         if mod(lev,2)==0 && lev < nlevel+1
%             
%             ids = abs(VZW-t(lev))<1e-8;
%             perturb = repmat([-1 1],1,floor(length(VZW(ids))/2));
%             if mod(length(VZW(ids)),2)==1
%                 perturb = [perturb -perturb(end)];
%             end
%             %             perturb = randn(size(VZW(ids)));
%             VZW(ids) = VZW(ids) + 2*h*perturb;
%         end
%     end
    
    %     VZW(ids) = VZW(ids) + h^2*randn(size(VZW(ids))); % for LSC - O(h^2 perturbation)
    
end

write_gmsh_file(filename,VXW,VYW,VZW,EToVW);




