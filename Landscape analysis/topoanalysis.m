function topo_info=topoanalysis(DEM, dx, sample_l, center, plots, fraction)
% Written by A Kubo 3/4/2020

%Inputs 
%   DEM     matrix of elevation data [row x col]
%   dx      x and y resolution, must be equal for now [m]
% (Optional arguments)
%   sample_area   length of side of area to sample over, must be less than 
%                 DEM sides, area=sample_l^2
%   center        [row col] of the center of the area 
%                 if the center yields an area that is outside the DEM area
%                 the script will shift it over 
%   plots         0 or 1. if 1 it will plot the hillshade image of the
%                 sampled DEM & the power spectrum
%  fraction       % of points to plot on graphs 


    [Ny, Nx]=size(DEM);
    DEM(isnan(DEM))=0;
    DEM(DEM<0)=0;
%% error checking     
    if nargin<3 
        sample_l=5000;
    end 
    
    if nargin<4
      center=[Ny/2, Nx/2];
    end 
    
    if nargin<5 
        plots=0 ;
    end 
    
    if nargin<6 
        fraction=0.1;
    end
    
    
    l= ceil(sample_l/dx);
    
    if Ny*dx < sample_l || Nx*dx < sample_l 
        msg='Error. Sample area is too large for DEM';
        error(msg)
    end 
    
    if center(2)<l
        center(2)=center(2)+abs((center(2)-l))+1;
    end 
    
    if center(2)+l>Nx 
        center(2)= center(2)- (center(2)+l - Nx);
    end 
    if center(1)+l>Ny 
        center(1)= center(1)- (center(1)+l - Ny);
    end 
    
    if center(1)<l
        center(1)=center(1)+abs((center(1)-l))+1;
    end 
    %% take a 5km x 5km section at the middle 
    l= ceil(sample_l/dx)/2;

    DEM_sample= DEM((center(1) -l):(center(1) +l), (center(2) -l):(center(2) +l));
    %new size
    [row, col]=size(DEM_sample);
    middle=[ceil(row/2), ceil(col/2)]
    %% detrend data 
    [DEM_detrend,~,~] = Detrend2(DEM_sample);

    % best fit plane
    plane=DEM_sample-DEM_detrend;

    % background slope
    dzdx = (plane(middle) - plane(middle(1)-1, middle(2)))/(dx);
    dzdy = (plane(middle) - plane(middle(1), middle(2)-1))/(dx);

    Slope_Gradient = sqrt(dzdx^2 + dzdy^2);
    %average slope in degrees
    gradient = atand(Slope_Gradient);
    %variance of detrended surface (m^2)
    variance = std(DEM_detrend(:))^2;

    pad=1;
    [Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);


    nbin=20;  % number of logarithmically spaced bins
    B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
    LowBin = 1;
    HighBin = nbin;

    % Fit trend to all bins
    fit = robustfit(B(:,1),B(:,2));
    beta = fit(2);

    %% how does it look 

    if plots
 
        h=hillshade(DEM, 1:Ny*dx, 1:Nx*dx, 'plotit');
        hold on 
        plot( [center(1)+l,center(1)+l], [center(2)-l,center(2)+l], 'w.-', 'LineWeight', 3)
        plot( [center(1)-l,center(1)-l], [center(2)-l,center(2)+l], 'w.-', 'LineWeight', 3)
        
        plot( [center(1)+l,center(1)-l], [center(2)-l,center(2)-l], 'w.-', 'LineWeight', 3)
        plot( [center(1)+l,center(1)-l], [center(2)+l,center(2)+l], 'w.-', 'LineWeight', 3)
        
        
        figure;
        %plot partial of spectrum
        RN = rand(1,length(Pfv));
        [~, index] = find(RN < fraction);
        P_frac = Pfv(index); 
        F_frac= ffv(index);
        loglog(F_frac, P_frac, '.' , 'Color', [0.5 0.5 0.5]);
        ylabel('DFT mean squared amplitude (m^2)')
        title('Periodogram')
        xlabel('Radial Frequency (1/m)')
        hold on 
        %plot fit 
        loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
        % plot binning 
        loglog(10.^B(LowBin:HighBin,1),10.^B(LowBin:HighBin,2),'ok','markerfacecolor','w', 'MarkerSize',9)
        
        %txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope)];
        %txt2= ["Varience = "  + num2str(ceil(varience))];
        
    end 
topo_info=[beta, variance, gradient];
end 
