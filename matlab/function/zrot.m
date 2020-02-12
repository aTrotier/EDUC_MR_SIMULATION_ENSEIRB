function Rz=zrot(phi,type)
%type : 'rad' or 'deg'
if nargin<2
    type='rad';
end

if type=='deg'
    Rz = [cosd(phi) -sind(phi) 0;sind(phi) cosd(phi) 0; 0 0 1];
else
    Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];
end
