function [npoin,nelem,nvfix,ntype,nnode,ndofn,ndime,ngaus, ...
	       nstre,nmats,nprop,lnods,matno,coord,props,nofix,  ...
	       iffix,fixed]=input_fem()

global in;
global out;


%read the input data:

npoin=fscanf(in,"%d",1);
nelem=fscanf(in,"%d",1);
nvfix=fscanf(in,"%d",1);
ntype=fscanf(in,"%d",1);
nnode=fscanf(in,"%d",1);
ndofn=fscanf(in,"%d",1); 
ndime=fscanf(in,"%d",1); 
ngaus=fscanf(in,"%d",1);
nstre=fscanf(in,"%d",1);
nmats=fscanf(in,"%d",1);
nprop=fscanf(in,"%d",1);


%read the element node numbers & material property number

for ielem=1:nelem
jelem=fscanf(in,"%d",1);
dummy=fscanf(in,"%d",[nnode+1,1]);
for inode=1:nnode
lnods(jelem,inode)=dummy(inode);
end
matno(jelem)=dummy(nnode+1);
end

%read nodal coordinates:

for ipoin=1:npoin
jpoin=fscanf(in,"%d",1);
dummy=fscanf(in,"%lf %lf",[2,1]);
for idime=1:ndime
coord(ipoin,idime)=dummy(idime);
end
end

%Read constraint nodes and their values:

for ivfix=1:nvfix
nofix(ivfix)=fscanf(in,"%d",1);
dummy1=fscanf(in,"%d %d",[2,1]);
dummy2=fscanf(in,"%lf %lf",[2,1]);
for idime=1:ndime
iffix(ivfix,idime)=dummy1(idime);
fixed(ivfix,idime)=dummy2(idime);
end
end

%read Material properties:
for imats=1:nmats
jmats=fscanf(in,"%d",1);
dummy=fscanf(in,"%lf %lf",[2,1]);
for iprop=1:nprop
props(jmats,iprop)=dummy(iprop);
end
end


%%printout:

fprintf(out,"************************ Results *****************************\n");

fprintf(out,"Number of Elements          : %5d\n",nelem);
fprintf(out,"Number of Node              : %5d\n",npoin);
fprintf(out,"Number of Fixed nodes       : %5d\n",nvfix);
fprintf(out,"Number of Nodes per element : %5d\n",nnode);
fprintf(out,"Number of Integration points: %5d\n",ngaus);
fprintf(out,"Number of Materials         : %5d\n",nmats);
fprintf(out,"Number of properties        : %5d\n",nprop);

if(ntype == 1)
  fprintf(out,"Plane-strain elasticity solution\n");
end
if(ntype == 2)
  fprintf(out,"Plane-stress elasticity solution\n");
end

fprintf(out,"Element connectivity\n");
fprintf(out,"Element No         Node numbers         Property Id\n");

for ielem=1:nelem

fprintf(out,"%5d",ielem);
for inode=1:nnode
fprintf(out,"%5d",lnods(ielem,inode));
end
fprintf(out,"%5d\n",matno(ielem));
end

fprintf(out,"Nodal Coordinates\n");
fprintf(out," Node Number     X-cord        Y-cord\n");
for ipoin=1:npoin
fprintf(out,"%5d   %14.6e  %14.6e\n", ipoin,coord(ipoin,1),coord(ipoin,2));
end

fprintf(out,"Fixed Nodes  Fixed DOFs    Values\n");

for ivfix=1:nvfix     

	     fprintf(out,"%5d %5d %5d %14.6e %14.6e\n",nofix(ivfix),iffix(ivfix,1),iffix(ivfix,2), ...
	fixed(ivfix,1),fixed(ivfix,2));
end

fprintf(out,"Material Properties\n")
for imats=1:nmats
fprintf(out," %d %14.6e %14.6e\n", props(imats,1),props(imats,2));
end





end %end function
