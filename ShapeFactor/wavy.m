
function im=wavy(r,a,f,name)
%wavy 
% written by A Kubo 4/2020
%% Inputs 
%   r - radius of circle in pixel units
%   a - amplitude of purtubation in pixel units
%   f - frequency of purtubation in pixel units
%       wavenlength will be 1/n

%% Outputs 
%   im - binary matrix of size resxres 
%   with 1s indicating the circle 

% creates wavey circles 
%dx=1/(4*f);
%   res - image resolution, output is resxres
%         res = 3*r/dx 
%         dx = 1/(4*f)
%res= 3*r/(dx);
res=100;
t=linspace(0,360,res);

x=(r+a*sind(f*t)).*(cosd(t));
y=(r+a*sind(f*t)).*(sind(t));

s=linspace(-1.5*r, 1.5*r, 3*r);
[xx,yy]=meshgrid(s,s);

im=inpolygon(xx,yy, x,y);

if nargin>3
    save(name, 'im', 'r', 'a', 'f')
end
%figure; imshow(im)