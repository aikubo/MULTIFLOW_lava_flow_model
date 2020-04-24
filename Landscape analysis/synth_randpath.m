clear all
close all
%% Generate Synthetic landscape
% ny=512;
% nx=512;
% dx=1;
% H=0.5;
% pad=0;
% periodic=1 ;
% var=600;
% beta=2+2*H;
% window=0;
% 
% [DEM Pm fm Pv fv] = synthspecNEW(nx,ny,dx,H,pad,window,var);
% a1=1 ;
% b1=0;
% c1=0;
%P = makeplane(nx,ny,a1,b1, c1 );
%DEM=M+P;

% steepest descent pathway
% [M, N] = size(DEM) ; % M x N : Y-dimension x X-dimension
% % fill DEM w/TopoToolbox
% DEMtopo = GRIDobj(1:N,1:M,DEM);
% DEMf = fillsinks(DEMtopo);
% % calculate drainage directions
% FD = FLOWobj(DEMf, 'preprocess');
% ST=STREAMobj(FD, flowacc(FD)>1000);
% S=trunk(ST);
%% steepest descent pathway
load stream_net.mat
% contains 
% Sm - STREAMobj 
% nx,ny - cols,rows 
% H - roughness parameter 
% var - varience 
% P -  plane 
% M - roughtopograph 
% DEM - M+P topography
% dx - spatial delta
DEMtopo = GRIDobj(1:ny,1:nx,DEM);


% % plot it
% figure
% subplot(1,2,1)
% imageschs(DEMtopo)%,[],'colormap',[1 1 1],'colorbar',false)
% hold on
% plot(Sm, 'r', 'LineWidth',3);
% subplot(1,2,2)
% plotdz(Sm,DEMtopo)



%% Curvature

% make rr for curvature code 
xx=Sm.x;
yy=Sm.y;

% use detrended or full elevation??
z=zeros(length(xx),1);
dist=zeros(length(xx),1);
z(1)=DEM(xx(1), yy(1));
for i=2:length(xx)
    z(i)=DEM(xx(i), yy(i));
    dist(i)=sqrt((xx(i)-xx(i-1)).^2 + (yy(i)-yy(i-1)).^2);
end 

rr=[xx'; yy' ;z'];

lwin=10;
weight=0;

[kappalist,taulist,ttlist,nnlist,bblist]=frenet_robust(rr,lwin,weight);
kappalist(isnan(kappalist))=0;



N=2^(nextpow2(length(rr)));
k = [0:N/2 (-N/2+1):-1] ;
f = k./(2*pi);
kappabar=fft(kappalist,N);

rrbar=fft(z,N);

% figure; 
% % subplot(2,1,1)
% % plot(cumsum(dist),y) 
% % xlabel('Distance (m)')
% % ylabel('Elevation (m)')
% % subplot(2,1,2)
% plot(cumsum(dist),kappalist, 'r')
% xlabel('Distance (m)')
% ylabel('Curvature (1/m)') 

figure; 
 subplot(1,2,1) 

P1=log((rrbar(1:N/2).*conj(rrbar(1:N/2))));

hold on 
nbin=30;

freq=log(f(1:N/2));
freq(freq<0)=0;
B1= bin(freq, P1, 30);
fit1=robustfit(B1(:,1), B1(:,8))
plot(freq, P1, 'k.')
hold on
plot(freq, log((f(1:N/2).^fit1(2)))+fit1(1), 'b', 'LineWidth', 3)
Spectral_slope1 = abs(fit1(2));
txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope1)];
text(1,20, txt, 'FontSize', 12)
xlabel('Frequency (1/m)')
ylabel('PSD')   
title('Topography transect')

subplot(1,2,2)

P2=log((kappabar(1:N/2).*conj(kappabar(1:N/2))));
B2= bin(freq, P2, 30);
fit2=robustfit(B2(:,1), B2(:,8))
Spectral_slope2 = abs(fit2(2));
txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope2)];

plot( freq, P2, 'k.')
hold on
plot(freq, log(f(1:N/2).^fit2(2))+fit2(1), 'r', 'LineWidth', 3)
text(1,6, txt, 'FontSize', 12)
xlabel('Frequency (1/m)')
ylabel('PSD')   
title('Curvature')

% figure; plot(k(1:N/2),abs(kappabar(1:N/2)))
% xlabel('wavenumber k')
% ylabel('PSD')
% 
