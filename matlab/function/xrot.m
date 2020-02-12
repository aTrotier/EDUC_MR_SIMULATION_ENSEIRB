function Rx=xrot(phi,type)
%type : 'rad' or 'deg'
if nargin<2
    type='rad';
end

if type=='deg'
    Rx = [1 0 0; 0 cosd(phi) -sind(phi);0 sind(phi) cosd(phi)];
else
    Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
end