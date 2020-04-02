
function im=wavy(r,a,n, res)
%wavy 
% written by A Kubo 4/2020
%% Inputs 
%   r - radius of circle in pixel units
%   a - amplitude of purtubation in pixel units
%   n - period of purtubation in pixel units
%   res - image resolution, output is resxres
%% Outputs 
%   im - binary matrix of size resxres 
%   with 1s indicating the circle 

% creates wavey circles 

t=linspace(0,360,res);

x=(r+a*sind(n*t)).*(cosd(t));
y=(r+a*sind(n*t)).*(sind(t));

s=linspace(-1.5*r, 1.5*r, res);
[xx,yy]=meshgrid(s,s);

im=inpolygon(xx,yy, x,y);
%figure; imshow(im)