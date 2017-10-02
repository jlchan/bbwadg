function minJ = min_pyr_jacobian(VX,VY,VZ)

% plotting points
[a b c] = meshgrid(linspace(-1,1,25));
a = a(:); b = b(:); c = c(:);
[rp sp tp] = pyr_abctorst(a,b,c);

[xp,yp,zp,~,~,~,~,~,~,~,~,~,Jp] = pyr_geom_factors(VX,VY,VZ,rp,sp,tp);

minJ = min(Jp);
if nargout==0
    color_line3(xp,yp,zp,Jp,'.')
    title(sprintf('Min wedge Jacobian = %f',minJ))
end

return
