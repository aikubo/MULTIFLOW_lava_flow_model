clear all
close all
%% Generate Synthetic landscape
ny=512;
nx=512;
dx=1;
H=0.5;
pad=0;
periodic=1 ;
var=300;
beta2d=2+2*H;
beta1d=1+2*H;
window=0;
H=linspace(0.1,1, 15);
C={};
beta=linspace(2,4,15);

for ii=1:length(H)
    [M1, ~, ~] = RedNoise(ny,nx, dx,beta(ii),var,periodic);
    a1=-1 ;
    b1=0;
    c1=500;
    P = makeplane(nx,ny,a1,b1, c1 );
    DEM=M1+P;

    [M, N] = size(DEM) ; % M x N : Y-dimension x X-dimension
    % fill DEM w/TopoToolbox
    DEMtopo = GRIDobj(1:N,1:M,DEM);
    DEMf = fillsinks(DEMtopo);
    % calculate drainage directions
    FD = FLOWobj(DEMf);
    ST=STREAMobj(FD, flowacc(FD)>1000);
    S=trunk(ST);

    Sm=modify(S, 'interactive', 'reachselect');
    close 
    
    C{ii,1}=H(ii);
    C{ii,3}=DEM
    C{ii,4}=Sm.x;
    C{ii,5}=Sm.y;
    
end

save('steepest_paths_var300.mat', 'C', 'P', 'ny', 'nx', 'dx', 'beta', 'var')
