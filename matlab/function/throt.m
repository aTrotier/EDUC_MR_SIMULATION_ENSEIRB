function Rth=throt(phi,theta,type)
%type : 'rad' or 'deg'
if nargin<3
    type='rad';
end

Rz = zrot(-theta,type);
Rx = xrot(phi,type);
Rth = inv(Rz)*Rx*Rz;


