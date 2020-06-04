function curve = naturalCurveD(k,t,isplotted)
%naturalCurveD reconstruct a space curve from curvature and torsion. 
%It uses Runge-Kutta method for solving ODE.
% The formula for updating tangent, normal and binormal vector defined in 
% a TNB frame can be found from reference:
%
% Sungyeop Lim and Soonhung Han, Helical extension method for solving the  
% natural equation of a space curve. 
% Surf. Topogr.: Metrol. Prop. 5 (2017) 035002
%
% However, the helical extension method introduced in the paper is not yet
% implemented. 
%
%   Input variable:
%   k(kappa)            - Curvature, which can be a single value or a
%                         vector
%   t(tau)              - Torsion, which can be a single value or a vector
%   isplotted           - binary value specifying whether the reconstructed
%                         curve is plotted or not.
%   Output variable                      
%   curve               - Reconstructed curve
%
% Author: Lei Xu, lei.xu@ki.se / lei.xu.optics@gmail.com, November 21, 2018


if nargin <3
    isplotted = 0;
end
ds = 2^(-nextpow2(max(k)))/10;% interval length
if numel(k)==1
    n = 100;
    k = k*ones(n,1);
    t = t*ones(n,1);
else
    n = length(k);
end
assert(n==length(t),'Curvature and torsion must have identical number of elements!');
T = zeros(n,3);
N = zeros(n,3);
B = zeros(n,3);
T(1,:) = [0 1 0];
B(1,:)=[1 0 0];
for i = 1:n
    if i >1
        T(i,:) = (1+k(i-1)^2*ds^2/2+(k(i-1)^4+k(i-1)^2*t(i-1)^2)*ds^4/4)*T(i-1,:)...
              + (k(i-1)*ds-(k(i-1)^3-k(i-1)*t(i-1)^2)*ds^3/6)*N(i-1,:)...
              + (k(i-1)*t(i-1)*ds^2/2-(k(i-1)^3*t(i-1)+k(i-1)*t(i-1)^3)*ds^4/24)*B(i-1,:);

        B(i,:) = (k(i-1)*t(i-1)*ds^2/2 - (k(i-1)*t(i-1)^3+k(i-1)^3*t(i-1))*ds^4/24)*T(i-1,:)...
              + (-t(i-1)*ds+(k(i-1)^2*t(i-1)+t(i-1)^3)*ds^3/6)*N(i-1,:)...
              + (1-t(i-1)^2*ds^2/2 + (k(i-1)^2*t(i-1)^2+t(i-1)^4)*ds^4/24)*B(i-1,:);

        N(i,:) = (-k(i-1)*ds + (k(i-1)*t(i-1)^2+k(i-1)^3)*ds^3/6)*T(i-1,:)...
              + (1-(k(i-1)^2+t(i-1)^2)*ds^2/2+(k(i-1)^2+t(i-1)^2)^2*ds^4/24)*N(i-1,:)...
              + (t(i-1)*ds - (k(i-1)^2*t(i-1)+t(i-1)^3)*ds^3/6)*B(i-1,:);  
    end  
    T(i,:) = T(i,:)/norm(T(i,:));
%     B(i,:) = B(i,:)/norm(B(i,:));
%     N(i,:) = cross(B(i,:),T(i,:));
    N(i,:) = cross(B(i,:),T(i,:))/norm(cross(B(i,:),T(i,:)));
    B(i,:) = cross(T(i,:),N(i,:));
    
end
% coordinates are calculated according to the trapezoidal rule.
curve = [0,0,0;cumsum((T(1:end-1,:)+T(2:end,:))/2*ds)];
if isplotted
    figure,
    plot3(curve(:,1),curve(:,2),curve(:,3))
    axis equal;
    axis vis3d,
    view([-50,10])
end
end

