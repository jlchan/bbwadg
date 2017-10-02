function minJ = min_tet_jacobian(VX,VY,VZ)

% plotting points
[a b c] = meshgrid(linspace(-1,1,25));
a = a(:); b = b(:); c = c(:);
[rp sp tp] = tet_abctorst(a,b,c);

[xp,yp,zp,~,~,~,~,~,~,~,~,~,Jp] = tet_geom_factors(VX,VY,VZ,rp,sp,tp);
if nargout ==0
    color_line3(xp,yp,zp,Jp,'.')
    title(sprintf('Min wedge Jacobian = %f',min(Jp)))
end
minJ = min(Jp);
return
