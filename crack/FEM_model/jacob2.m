function [cartd,djacb,gpcod] =jacob2(ielem,elcod,kgasp,shape,deriv,nnode,ndime)

format long;

%gauss point coordinates:

cg = elcod*shape';

gpcod(1,kgasp)=cg(1);
gpcod(2,kgasp)=cg(2);

%jacobian

xjacm = deriv*elcod';

%Determinate of Jacobian

djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1);

if(djacb <= 0.0)
fprintf('Element No: %5d\n',ielem);
error('Program terminated zero or negative area');
end

%cartesion derivatives:

xjaci(1,1)=xjacm(2,2)/djacb;
xjaci(2,2)=xjacm(1,1)/djacb;
xjaci(1,2)=-xjacm(1,2)/djacb;
xjaci(2,1)=-xjacm(2,1)/djacb;

cartd = xjaci*deriv;


end %endfunction
