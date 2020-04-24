clear all
close all
%% Generate Synthetic landscape
load steepest_paths_var300.mat

nruns=15;
curv=zeros(nruns,1);
descentx=zeros(nruns,1);
descenty=zeros(nruns,1);
transect=zeros(nruns,1);
topo2D=zeros(nruns,1);
beta=zeros(nruns,1);

for i=1:nruns 
    beta(i)=C{i,1};
end 


for ii=1:nruns
    DEM=C{ii,3};
    xx=C{ii,4};
    yy=C{ii,5};
    
    %% 2d spectral 
    DEM_detrend=DEM-P; 
    pad = 0; % 1 means pad the data with zeros to a power of 2, 0 no padding
    window = 1; % 1 means window the data prior to taking the FFT, 0 no window

    [Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);
    nbin=30;
    B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
    LowBin = 1;
    HighBin = nbin;

    % Fit trend to all bins
    fit2D = robustfit(B(:,1),B(:,2));
    Spectral_slope_2D = fit2D(2);

    %% TOPO2D FITS ANALYTICAL 4/23
    
    %% pathways
    % use detrended or full elevation??
    z=zeros(length(xx),1);
    z_detrend=zeros(length(xx),1);
    dist=zeros(length(xx),1);
    z(1)=DEM(xx(1), yy(1))-P(xx(1), yy(1)) ;
    for i=2:length(xx)
        z(i)=DEM(xx(i), yy(i));
        z_detrend(i)=DEM(xx(i), yy(i))-P(xx(i), yy(i));
        dist(i)=sqrt((xx(i)-xx(i-1)).^2 + (yy(i)-yy(i-1)).^2 + (z(i)-z(i-1)).^2);
    end 

    %% curvature
    rr=[xx'; yy' ;z'];

    lwin=10;
    weight=0;

    [kappalist,taulist,ttlist,nnlist,bblist]=frenet_robust(rr,lwin,weight);
    kappalist(isnan(kappalist))=0;


%% fft of profiles

%% curvature
    N=2^(nextpow2(length(rr)));
    k = [0:N/2 (-N/2+1):-1] ;
    f = k./(2*pi);
    kappabar=fft(kappalist,N);
    
    Pkappa=log((kappabar(1:N/2).*conj(kappabar(1:N/2))));

    nbin=30;

    freq=log(f(1:N/2));
    freq(freq<0)=0;
    B1= bin(freq, Pkappa, 30);
    fit1=robustfit(B1(:,1), B1(:,8));
    Spectral_slopekappa = abs(fit1(2));
    

%% 1D profile
    col=nx/2;
    nprofiles=100;
    slopeprofile=0;
    N=nx;
    for bb=1:nprofiles
    
        z_transect= DEM_detrend(:,col+bb);
        transectbar=fft(z_transect,N);
        Ptrans=log((transectbar(1:N/2).*conj(transectbar(1:N/2)))/N);
        k = [0:N/2 (-N/2+1):-1] ;
        f = k./(2*pi);
        freq=log(f(1:N/2));
        freq(freq<0)=0;
        B2= bin(freq, Ptrans, 30);
        fit2=robustfit(B2(:,1), B2(:,8));
        Spectral = abs(fit2(2));
        slopeprofile = slopeprofile+Spectral;
    end 
    
    Spectral_slopeprofile=slopeprofile/nprofiles;
    %% TOPO1D FITS ANALYTICAL 4/23
    
%% steepest descent
    n1=length(xx);
    N=2^nextpow2(n1);
    xx=xx(:);
    pathxbar=fft(xx,N);
    xbar=log(pathxbar(1:N/2).*conj(pathxbar(1:N/2)));
    k = [0:N/2 (-N/2+1):-1] ;
    f = k./(2*pi);
    freqt=log(f(1:N/2));
    freqt(freqt<0)=0;
    Bx= bin(freqt, xbar, 30);
    yy=yy(:);
    
    
    pathybar=fft(yy,N);
    
    ybar=log(pathybar(1:N/2).*conj(pathybar(1:N/2)));

    k = [0:N/2 (-N/2+1):-1] ;
    f = k./(2*pi);
    freqt=log(f(1:N/2));
    freqt(freqt<0)=0;
    Bx= bin(freqt, xbar, 30);
    By= bin(freqt, ybar, 30);
    fitx=robustfit(Bx(:,1), Bx(:,8));
    fity=robustfit(By(:,1), By(:,8));
    Spectral_slopex = abs(fitx(2));
    Spectral_slopey = abs(fity(2));
    
    
    % figure; 
    % % subplot(2,1,1)
    % % plot(cumsum(dist),y) 
    % % xlabel('Distance (m)')
    % % ylabel('Elevation (m)')
    % % subplot(2,1,2)
    % plot(cumsum(dist),kappalist, 'r')
    % xlabel('Distance (m)')
    % ylabel('Curvature (1/m)') 



%     figure; 
%      subplot(1,2,1)     
%     plot(freq, P1, 'k.')
%     hold on
%     plot(freq, log((f(1:N/2).^fit1(2)))+fit1(1), 'b', 'LineWidth', 3)
%     txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope1)];
%     text(1,20, txt, 'FontSize', 12)
%     xlabel('Frequency (1/m)')
%     ylabel('PSD')   
%     title('Topography transect')



%     subplot(1,2,2)    
%     plot( freq, P2, 'k.')
%    txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope2)];
%     hold on
%     plot(freq, log(f(1:N/2).^fit2(2))+fit2(1), 'r', 'LineWidth', 3)
%     text(1,6, txt, 'FontSize', 12)
%     xlabel('Frequency (1/m)')
%     ylabel('PSD')   
%     title('Curvature')


    


    descentx(ii)=Spectral_slopex;
    descenty(ii)=Spectral_slopey;
    curv(ii)=Spectral_slopekappa;
    transect(ii)=Spectral_slopeprofile;
    topo2D(ii)=Spectral_slope_2D;
    % figure; plot(k(1:N/2),abs(kappabar(1:N/2)))
    % xlabel('wavenumber k')
    % ylabel('PSD')
    % 
end

save('fft_topo_profile_descent_curv.mat', 'beta','topo2D','transect', 'descentx', 'descenty','curv')
