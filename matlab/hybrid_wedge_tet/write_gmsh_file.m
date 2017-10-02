% writes gmsh file given VXYZ and EToV. works with TET and WEDGE only!

function write_gmsh_file(filename,VX,VY,VZ,EToV)

gmshfile = ['Grid/' filename '.msh'];
fprintf(sprintf(gmshfile))
fid = fopen(sprintf(gmshfile), 'w');
fprintf(fid, '$MeshFormat\n');
fprintf(fid, '2.2 0 8\n');
fprintf(fid, '$EndMeshFormat\n');
% 
fprintf(fid, '$Nodes\n');
fprintf(fid, '%d\n',length(VX));
for i = 1:length(VX)
    fprintf(fid,'%d %f %f %f\n',i,VX(i),VY(i),VZ(i));
end
fprintf(fid,'$EndNodes\n');
% 
fprintf(fid, '$Elements\n');
fprintf(fid, '%d\n',size(EToV,1));
for e = 1:size(EToV,1)
    v = EToV(e,:); 
    nverts = nnz(v>0);
    if nverts==4 % tet
        type = 4; 
    elseif nverts==6 % wedge
        type = 6;
    else % uhh....
        disp('Error: wrong # of vertices here?  Only tet/wedges supported.')
    end
    if type==6
        fprintf(fid,'%d %d 2 200 20 ',e,type); % elem #, elem type, 0 tags
    elseif type==4
        fprintf(fid,'%d %d 2 400 20 ',e,type); % elem #, elem type, 0 tags
    end
    for v = 1:nverts
        fprintf(fid,'%d ',EToV(e,v));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'$EndElements\n');
fclose(fid)

disp(['wrote gmsh file ' sprintf(gmshfile)])