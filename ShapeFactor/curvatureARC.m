FlowMap=Kilauea_shape;
plotson=1; 
close all
%function [rr_con, kappalist_grid, kappalist_con, B, Spectral_Slope, intercept]=curvatureARC(FlowMap, VentLocation, plotson)
[row,col]=size(FlowMap);

[xs,ys]=meshgrid(1:col, 1:row);


%% Find perimeter 
EDGES=bwboundaries(FlowMap);
Xloc=EDGES{1}(:,1); % ordered set of boundary points X, columns
Yloc=EDGES{1}(:,2); % ordered set of boundary points y, rows

%Xloc=xs(Xloc);
%Yloc=ys(Yloc);

%% Make Vent the starting point 
VentX=VentLocation(1);
VentY=VentLocation(2);
boundary=length(Xloc);
dvx= VentX-Xloc;
dvy=VentY-Yloc;
Vdist=sqrt( (dvx).^2 + (dvy).^2);
vboundary=find(Vdist==min(Vdist));

%make sure that there is only one closest point 
if length(vboundary)>1
    vboundary=vboundary(1);
end

% shift so that in the 1st position is the point closest to the vent
% now poinst are in a ordered perimeter (Xloc,Yloc) starting at the vent
Xloc=circshift(Xloc, (boundary-vboundary));
Yloc=circshift(Yloc, (boundary-vboundary));



dlist=[];
dlist(1)=0; 
for i=length(Xloc):-1:2
    dlist(i)=sqrt( (Xloc(i)-Xloc(i-1))^2 + (Yloc(i)-Yloc(i-1))^2);
end 

z=zeros(1,length(Xloc));
rr_grid=[z; Xloc(:)'; Yloc(:)'];


%% calculate the Frenet-Serret frame 
%now we can calculate curvature
weight=0.0; %regularization weight
lwin=10; %moving window defining the spatial resolution of curvature
% compute the Frenet frame ===============================================
[kappalist_grid,taulist,ttlist,nnlist,bblist]=frenet_robust(rr_grid,lwin,weight);
kappalist_grid(isnan(kappalist_grid))=0;
taulist(isnan(taulist))=0;
ttlist(isnan(ttlist))=0;
bblist(isnan(bblist))=0;
nnlist(isnan(nnlist))=0;
% figure;
% plot(cumsum(dlist), ttlist(2,:), '*')
% hold on 
ds=mean(dlist);

[tt2,d]=resample(ttlist(2,:),cumsum(dlist), 1/ds);
[tt3,~]=resample(ttlist(3,:),cumsum(dlist), 1/ds);

%

rr_con=zeros(3,length(tt2));
rr_con(:,1)= rr_grid(:,1);
for i=2:length(tt2)
    rr_con(2,i)= rr_con(2,i-1)+(tt2(i)*ds);
    rr_con(3,i)= rr_con(3,i-1)+(tt3(i)*ds);
    %rr_con(3,i)= rr_con(3,i-1)+(ttlist(i)*ds);
end

[kappalist_con,~,~,~,~]=frenet_robust(rr_con,lwin,weight);
L = sum(dlist);

kappalist_con(isnan(kappalist_con))=0;

N_grid=2^nextpow2(length(rr_grid));
kappa_grid=fft(kappalist_grid,N_grid);
k_grid = [0:N_grid/2 (-N_grid/2+1):-1] *2*pi /L;

N_con=2^nextpow2(length(rr_con));
kappabar_con=fft(kappalist_con,N_con);
k_con = [0:N_con/2 (-N_con/2+1):-1] *2*pi /L;
nbin=20;
f=log10(k_con(1:N_con/2));
P=log10(abs(kappabar_con(1:N_con/2)));
P(isnan(P))=0;
f(isnan(f))=0;
f(isinf(f))=0;
B = bin( f, P,nbin,0); % bin the log-transformed data. 
bcol=8;
% % Fit trend to all bins
fit = robustfit(B(12:18,1),B(12:18,8));
Spectral_Slope = fit(2);
intercept=fit(1);

if plotson
    figure;
    nfig=3;
    subplot(1,nfig,1)

    hold on 
    plot(rr_grid(2,:), rr_grid(3,:), 'r')
    plot(rr_con(2,:), rr_con(3,:), 'b')
    legend('Gridded', 'Reconstructed')
    axis equal

    set(gca,'fontsize',18)
    subplot(1,nfig,2)
    hold on 
    plot(cumsum(dlist), kappalist_grid, 'r')
    plot(d,kappalist_con, 'b', 'LineWidth',2)
    legend('Gridded', 'Reconstructed')
    set(gca,'fontsize',18)
    subplot(1,nfig,3)
    hold on 
    plot( log10(k_grid(1:N_grid/2)), log10(abs(kappa_grid(1:N_grid/2))), '.', 'Color',[1 .6 .6])
    plot(log10(k_con(1:N_con/2)), log10(abs(kappabar_con(1:N_con/2))),'b.', 'LineWidth',2)
    
    plot( B(:,1), (intercept)+ Spectral_Slope*B(:,1), 'k')
    plot(B(11:19,1), B(11:19,8), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black')
    
    legend('Gridded', 'Reconstructed')
    set(gca,'fontsize',18)
    
end 

