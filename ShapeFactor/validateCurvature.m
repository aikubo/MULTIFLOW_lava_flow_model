load wavy_10.mat
dx=1;
VentLocation =[75,75];
unwindFlow=curvatureNEW(im, VentLocation, 0,0);
%plotCurvature(im, VentLocation, unwindFlow);
Perimeter=max(unwindFlow(:,1));
L=length(unwindFlow);
%% Analytical Solution 
t=linspace(0,360,L);

%x=(r+a*sind(f*t)).*(cosd(t));
%y=(r+a*sind(f*t)).*(sind(t));
r_a=r+a*sind(f.*t);


%%
% frequency 
% i set this already 
fUnwind=f/Perimeter;
Fs=L/Perimeter;
L=length(unwindFlow);
n = 2^nextpow2(L);

X=fft(unwindFlow(:,1), n);


% P2=abs(FTFlow/L);
% P1=P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
x= unwindFlow(:,3);
np=3;
figure;
subplot(1,np,1) 
imshow(im)
hold on 
plot(75,75, 'r+', 'LineWidth', 2)
subplot(1,np,2)
plot(t, unwindFlow(:,3), 'LineWidth', 1);
hold on 
plot(t,r_a)
ylim([35,65])
legend( 'UnwindFlow', 'Analytical')

title('Space Domain')

% subplot(1,np,2)
% NFFT=1024; %NFFT-point DFT      
% X=fft(x,NFFT); %compute DFT using FFT        
% nVals=0:NFFT-1; %DFT Sample points       
% plot(nVals,abs(X));      
% title('Double Sided FFT - without FFTShift');        
% xlabel('Sample points (N-point DFT)')        
% ylabel('DFT Values');
% 
% subplot(1,np,3)
% nVals=(0:NFFT-1)/NFFT; %Normalized DFT Sample points
% plot(nVals,abs(X));
% title('Double Sided FFT - without FFTShift');
% xlabel('Normalized Frequency')
% ylabel('DFT Values');
% 
% subplot(1,np,4)
% NFFT=1024; %NFFT-point DFT
% X=fftshift(fft(x,NFFT)); %compute DFT using FFT
% fVals=(-NFFT/2:NFFT/2-1)/NFFT; %DFT Sample points
% plot(fVals,abs(X));
% title('Double Sided FFT - with FFTShift');
% xlabel('Normalized Frequency')
% ylabel('DFT Values');
% 
% subplot(1,np,5)
% 
% fVals=Fs*(-NFFT/2:NFFT/2-1)/NFFT;
% plot(fVals,abs(X),'b');
% title('Double Sided FFT - with FFTShift');
% xlabel('Frequency (Hz)')
% ylabel('|DFT Values|');

subplot(1,np,3)
L=length(x);        
NFFT=1024;       
X=fft(x,NFFT);  

X2=fft(r_a, NFFT);

Px=X.*conj(X)/(NFFT*L); %Power of each freq components       
fVals=Fs*(0:NFFT/2-1)/NFFT;     

Px2=X2.*conj(X2)/(NFFT*L); %Power of each freq components       
fVals=Fs*(0:NFFT/2-1)/NFFT;      
plot(fVals,Px(1:NFFT/2),'b','LineWidth',1);   
hold on 
plot(fVals,Px2(1:NFFT/2),'r.-','LineWidth',1);    
xline(fUnwind, 'r','LineWidth',2); 
xlim([0.01, 0.03])
title('One Sided Power Spectral Density');       
xlabel('Frequency (Hz)')         
ylabel('PSD');
