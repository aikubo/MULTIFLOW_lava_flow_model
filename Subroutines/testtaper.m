%% testing tapering etc of synthetics 

load steepest_paths_var100.mat

ny=512;
nx=512;
dx=1;
pad=0;
periodic=1 ;
var=600;
nruns=15;

for ii=1:nruns
    H(ii)=C{ii,1};
end 
beta2d=2+2*H;
beta1d=1+2*H;
window=0;

topo2D_perron=zeros(nruns,1);
topo2D_kubo=zeros(nruns,1);
H=zeros(nruns,1);

for ii=1:nruns
    H(ii)=C{ii,1};
    DEM=C{ii,3};

    DEM_detrend=DEM-M1; 
        %% 2d PERON spectral 
    pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
    window = 1; % 1 means window the data prior to taking the FFT, 0 no window

    [Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);
    nbin=30;
    B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
    LowBin = 1;
    HighBin = nbin;

    % Fit trend to all bins
    fit2D = robustfit(B(:,1),B(:,2));
    Spectral_slope_2D = fit2D(2);
    
    topo2D_perron(ii)=Spectral_slope_2D;
    
    %% 2D Kubo spectral 
    
    [ny,nx]=size(DEM_detrend);
    Ly=ny;
    Lx=nx;
    dfx = 1/(dx*Lx);
    dfy = 1/(dx*Ly);

    M=fftshift(fft2(fftshift(DEM_detrend), ny, nx));
    M(Ly/2 + 1, Lx/2 + 1)=0;
    Wss = sum(sum(ones([ny nx]))); 
    M= M .* conj(M); %/(ny*nx*Wss);
    Pm=M;
    % Create a matrix of radial frequencies
    xc = Lx/2+1; yc = Ly/2+1; % matrix indices of zero frequency
    [cols rows] = meshgrid(1:Lx,1:Ly); % matrices of column and row indices
    fm = sqrt((dfy*(rows-yc)).^2 + (dfx*(cols-xc)).^2); % frequency matrix

    % Create sorted, non-redundant vectors of frequency and power 
    M = M(:,1:(Lx/2+1));
    % fvec = dfreq*sqrt((rows(:,1:xc)-yc).^2 + (cols(:,1:xc)-xc).^2);
    fv = fm(:,1:(Lx/2+1));

    fv((yc+1):Ly,xc) = -1; % This half-column is redundant. Set the 
                           % frequencies to negative values so they 
                           % will be clipped out below
    fv = sortrows(horzcat(fv(:),M(:)),1); % concatenate column vectors of 
                                          % frequency and PSD and sort by 
                                          % frequency
    fv = fv(fv(:,1)>0,:); % Take only positive frequencies. This gets rid 
                          % of DC (zero freq) as well as the redundant 
                          % frequencies we set to -1 above

    % Separate into power and frequency vectors and assign to output arguments
    Pv = 2*fv(:,2); % the factor of 2 corrects for the fact that we have
                    % taken only half of the 2D spectrum. sum(Pvec) should
                    % now equal sum(Pmat(:)).
    fv = fv(:,1);
    
    B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 
    LowBin = 1;
    HighBin = nbin;

    % Fit trend to all bins
    fit2D = robustfit(B(:,1),B(:,2));
    Spectral_slope_2D = fit2D(2);
    
    topo2D_kubo(ii)=Spectral_slope_2D;
end


beta1D=1+2.*H;
beta2D=2+2.*H;

figure; 
plot(H,H, 'k')
hold on
plot(H, (abs(topo2D_perron)-2)/2, 'g.', 'MarkerSize', 30)
plot(H, (abs(topo2D_kubo)-2)/2, 'r+', 'MarkerSize', 30)

xlabel('Roughness Parameter')
ylabel('Spectral Slope of 2D topography')
legend('Analytical', 'Synthetic', 'Without Taper')
