function [dmatx] = modps(mtype,ntype,nstre,props)

format long;

%--- Material Parameters:

young=props(mtype,1);
poiss=props(mtype,2);

for istre=1:3
for jstre=1:3
dmatx(istre,jstre)=0.0;
end
end

if(ntype == 1)

%--- Plane Stress:

const=young/(1.0-poiss*poiss);
dmatx(1,1)=const;
dmatx(2,2)=const;
dmatx(1,2)=const*poiss;
dmatx(2,1)=const*poiss;
dmatx(3,3)=(1.0-2.0*poiss)*const/2.0;
end

if(ntype == 2)

%--- Plane Strain:

const=young*(1.0-poiss)/((1+poiss)*(1.0-2.0*poiss));
dmatx(1,1)=const;
dmatx(2,2)=const;
dmatx(1,2)=const*poiss/(1.0-poiss);
dmatx(2,1)=const*poiss/(1.0-poiss);
dmatx(3,3)=(1.0-2.0*poiss)*const/(2.0*(1.0-poiss));

end
end %endfunction

