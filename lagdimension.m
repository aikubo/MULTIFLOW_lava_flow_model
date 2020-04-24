load steepest_paths_var100.mat

ny=512;
nx=512;
dx=1;
pad=0;
periodic=1 ;
var=600;
nruns=15;
H=zeros(nruns,1);
for ii=1:nruns
    H(ii)=C{ii,1};
end 
beta2d=2+2*H;
beta1d=1+2*H;
D2=3-H;
H_calc=zeros(1,nruns);
Dlag=zeros(1,nruns);
window=0;

for ii=1:nruns

    DEM=C{ii,3};
    lag=[10:1:256];
    E=zeros(1,length(lag));
    
    [xx, yy]=meshgrid(1:nx, 1:ny);
    E(1)=0;

    for i=2:length(lag)
        xpt=randi([1 512]);
        ypt=randi([1 512]);
        dist= sqrt( (xx-xpt).^2 + (yy-ypt).^2); 
        dist=dist-lag(i);
        ind=find(dist==min(dist));
        pt=randi([1,length(ind)]);
        xpt2=xx(ind(pt));
        ypt2=yy(ind(pt));
        E(i)=(DEM(xpt2,ypt2)-DEM(xpt,ypt))^2+E(i-1);
    end 


    nbin=30;
    B = bin(log10(lag),log10(E),nbin,0); % bin the log-transformed data. 

    startbin=5;
    endbin=nbin-5;
    % Fit trend to all bins
    fit2D = robustfit(B(startbin:endbin,1),B(startbin:endbin,8));
    
    lagfit = fit2D(2);
    C_fit=fit2D(1);
    
    H_calc(ii)=(lagfit/2);
    Dlag(ii)=3-H_calc(ii);

    
    figure; 
    plot(log10(lag), log10(E), 'k.')
    hold on
    plot(B(startbin:endbin,1),B(startbin:endbin,8), 'r.', 'MarkerSize', 20)
    plot(log10(lag), lagfit*log10(lag)+log10(C_fit), 'g-')
end
beta_calc=7-2*Dlag;

figure;
plot(H,beta2D, 'k')
hold on 
plot(H, beta_calc, 'r+')