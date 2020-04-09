f=10; %frequency of sine wave
overSampRate=30; %oversampling rate
fs=overSampRate*f; %sampling frequency
phase = 1/3*pi; %desired phase shift in radians
nCyl = 5; %to generate five cycles of sine wave
 
t=0:1/fs:nCyl*1/f; %time base
np=6;  
x=sin(2*pi*f*t+phase); %replace with cos if a cosine wave is desired
figure; 
subplot(1,np,1)
plot(t,x);
title(['Sine Wave f=', num2str(f), 'Hz']);
xlabel('Time(s)');
ylabel('Amplitude');

subplot(1,np,2)
NFFT=1024; %NFFT-point DFT      
X=fft(x,NFFT); %compute DFT using FFT        
nVals=0:NFFT-1; %DFT Sample points       
plot(nVals,abs(X));      
title('Double Sided FFT - without FFTShift');        
xlabel('Sample points (N-point DFT)')        
ylabel('DFT Values');

subplot(1,np,3)
nVals=(0:NFFT-1)/NFFT; %Normalized DFT Sample points
plot(nVals,abs(X));
title('Double Sided FFT - without FFTShift');
xlabel('Normalized Frequency')
ylabel('DFT Values');

subplot(1,np,4)
NFFT=1024; %NFFT-point DFT
X=fftshift(fft(x,NFFT)); %compute DFT using FFT
fVals=(-NFFT/2:NFFT/2-1)/NFFT; %DFT Sample points
plot(fVals,abs(X));
title('Double Sided FFT - with FFTShift');
xlabel('Normalized Frequency')
ylabel('DFT Values');

subplot(1,np,5)

fVals=fs*(-NFFT/2:NFFT/2-1)/NFFT;
plot(fVals,abs(X),'b');
title('Double Sided FFT - with FFTShift');
xlabel('Frequency (Hz)')
ylabel('|DFT Values|');

subplot(1,np,6)
L=length(x);        
NFFT=1024;       
X=fft(x,NFFT);       
Px=X.*conj(X)/(NFFT*L); %Power of each freq components       
fVals=fs*(0:NFFT/2-1)/NFFT;      
plot(fVals,Px(1:NFFT/2),'b','LineSmoothing','on','LineWidth',1);         
title('One Sided Power Spectral Density');       
xlabel('Frequency (Hz)')         
ylabel('PSD');