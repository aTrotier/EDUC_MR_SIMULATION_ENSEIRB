function Ry=yrot(phi,type)
%type : 'rad' or 'deg'
if nargin<2
    type='rad';
end

if type=='deg'
    Ry = [cosd(phi) 0 sind(phi);0 1 0;-sind(phi) 0 cosd(phi)];
else
    Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];
end