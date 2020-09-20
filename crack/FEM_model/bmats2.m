function [bmatx]=bmats2(dgdx,nelem,nnode,nstre,nevab,kgasp)

format long;

bmatx = zeros(nelem,nstre,nevab);

ngash=0;

for inode=1:nnode

mgash=ngash+1;
ngash=mgash+1;

bmatx( :,1,mgash)=dgdx( :,kgasp,1,inode);
bmatx( :,1,ngash)=0.0;
bmatx( :,2,mgash)=0.0;
bmatx( :,2,ngash)=dgdx( :,kgasp,2,inode);
bmatx( :,3,mgash)=dgdx( :,kgasp,2,inode);
bmatx( :,3,ngash)=dgdx( :,kgasp,1,inode);

end
end %endfunction
