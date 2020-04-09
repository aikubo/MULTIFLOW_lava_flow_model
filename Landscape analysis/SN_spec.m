%Spectral Characteristics of Sierra NEgra for 
% 5km x 5km box near Sylvia's flow

%first window slightly to make detrending by a plane less bad...
SNDEM_spec=SNDEM(270:699,100:509); 
ysnS=ysnX(270:699); xsnS=xsnX(100:509);

%distance vectors in m
Xm=0:p.dx:p.dx*(length(xsnS)-1);Ym=0:p.dy:p.dy*(length(ysnS)-1);
%detrend
Zde=Detrend2(SNDEM_spec);
%background
ZBACK = (SNDEM_spec - Zde);

dzdx = (ZBACK(end,end) - ZBACK(end,end-1))/p.dx;
dzdy = (ZBACK(end,end) - ZBACK(end-1,end))/p.dy;

Slope_Gradient = sqrt(dzdx^2 + dzdy^2);
%average slope in degrees
Slope_Degree = atand(Slope_Gradient);
%variance of detrended surface (m^2)
Variance = std(Zde(:))^2;

[p.Ny, p.Nx] = size(SNDEM_spec);

%now compute spectra - doing this on windowed and tapered but NOT detrended DEM 
%this may be something to play with in the future
DEMtp = DEMtaper(Zde,p.Ny,p.Nx);

[Pmat, fmat, Pvec, fvec] = fft2D(Zde,p.dx,p.dx,1,0); 

%plot only a fraction of the DEM points for efficiency
FractionPlot = 0.1; %0.01; 
RN = rand(1,length(Pvec));
[~, index] = find(RN < FractionPlot);
PvecPLOT = Pvec(index); 
fvecPLOT = fvec(index);

% Plot the raw and binned versions of the 1D spectrum
f_PowerSpec = figure; 
hold on
% raw data
plot(fvecPLOT,PvecPLOT,'.','color',[0.5 0.5 0.5])

% Create a binned "1D" radial power spectrum
nbin = 12;  % number of logarithmically spaced bins
B = bin(log10(fvec),log10(Pvec),nbin,0); % bin the log-transformed data. 

plot(10.^B(:,1),10.^B(:,2),'ok','markersize',10, 'markerfacecolor','w')

% Fit trend to all bins
LowBin = 1;
HighBin = 12;

fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_all = fit(2);

disp(['Spectral slope in vicinity of Sylias flow is ' num2str(Spectral_slope_all)])
    
set(gcf,'color','w');
set(gca,'fontname','Times')
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'LineWidth',2);

set(f_PowerSpec, 'units', 'inches', 'pos', [0 0.5 7 3.5])
set(gca,'fontsize',14);
ylabel('Spectral power (m^2)','fontsize',18)
xlabel('Frequency (1/m)','fontsize',18)
title('S.N. TanDemX topographic power spectrum in 5 x 5 km box')
%clean up workspace
%clear Pvec fvec PvecPLOT fvecPLOT RN ZBACK Zde fmat Pmat
end
%% ---------------------------- DETREND and TAPER DEM -------------------------------
%not detrending for now
% DEM_detrend = Detrend2(SNDEM);
% DEM_plane = SNDEM - DEM_detrend; 
[p.Ny, p.Nx] = size(SNDEM);

DEMtp = DEMtaper(SNDEM,p.Ny,p.Nx);

%% ------------------------ LOW PASS FILTER DEM ---------------------------

% Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -      
% flo < fhi 
flo = 1/(FilteredWavelength + p.dx); % can modify as required   
fhi = 1/(FilteredWavelength); % can modify as required


% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMtp = SpecFilt2D(DEMtp,p.dx,p.dy,[flo fhi],'lowpass');

%trim to original size
SNDEM = DEMtp(p.Ny/2+1:3*p.Ny/2,p.Nx/2+1:3*p.Nx/2);
% Add the best-fit plane to the lowpass filtered landscapes
%SNDEM = SNDEM + DEM_plane;
SNDEM(isnan(SNDEM)==1) = nan; 

