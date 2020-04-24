%% Rednoise spectrum generation 1D 
% Written by AKH April 2020 
%
% Generation of 1D red noise spectrum 

% beta - slope of the rednoise 
% N    - number of points 
close all 

beta=2;
N=512; 
dx=5 % same as sampling freqency

% generate random ( white noise) 
x=randn(1,N);

% take ft
Xbar=fft(x);

% P=abs(Xbar).*2;
% freq=1/dx*(0:(N/2))/N;
% PSD=P(1:(N/2)+1); 
% plots=0;
% 
% if plots
%     figure;
%     subplot(1,2,1)
%     plot(t, x)
%     subplot(1,2,2)
%     plot(freq, PSD)
% end 

L=N/2+1; % unique fft points

f=(1:L)*1/dx; % frequency 
X=X(1:L); 
X=X./f;

X=[X conj(X(end-1:-1:2))];

plot(log10(f),log10(X(1:L)))

y=real(ifft(X));

