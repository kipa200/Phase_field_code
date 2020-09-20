function [ ] = write_vtk_fem(npoin,nelem,nnode,lnods,coord,istep,cont1)

format long;

%-- open file

fname=sprintf('time_%d.vtk',istep);
out =fopen(fname,'w');


%--- start writing:

%% header
fprintf(out,'# vtk DataFile Version 2.0\n');
fprintf(out,'time_10.vtk\n');
fprintf(out,'ASCII\n');
fprintf(out,'DATASET UNSTRUCTURED_GRID\n');

%write nodal coordinates:

fprintf(out,'POINTS %5d float\n',npoin);

dummy = 0.0;

for ipoin=1:npoin
fprintf(out,'%14.6f %14.6f %14.6f\n',coord(ipoin,1),coord(ipoin,2), dummy);
end

%--- write element connectivity:

iconst1 = nelem*(nnode+1);

fprintf(out,'CELLS %5d %5d\n', nelem,iconst1);

for ielem=1:nelem
fprintf(out,'%5d',nnode);
for inode=1:nnode
fprintf(out,'%5d',(lnods(ielem,inode)-1));
end
fprintf(out,'\n');
end

%--- write cell types:

if(nnode == 8)
ntype = 23;
end

if(nnode == 4)
ntype = 9;
end

if(nnode == 3)
ntype = 5;
end

fprintf(out,'CELL_TYPES %5d\n',nelem);

for i=1:nelem
fprintf(out,'%2d\n', ntype);
end


%--- write Nodal scaler & vector values:

fprintf(out,'POINT_DATA %5d\n',npoin);

%--- write concentration values as scalar:

fprintf(out,'SCALARS Con  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

for ipoin=1:npoin
fprintf(out,'%14.6e\n',cont1(ipoin));
end


fclose(out);

end %endfunction
