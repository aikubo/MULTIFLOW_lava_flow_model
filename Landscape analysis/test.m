close all

beta=3.2
M=512;
% signal parameters
fs = 44100;         % sampling frequency, Hz
T = 60;             % signal duration, s
N = round(fs*T);    % number of samples
% noise generation
M=N;
% generate white noise sequence
x = randn(1, M);
% FFT
X = fft(x);
% prepare a vector with frequency indexes 
NumUniquePts = M/2 + 1;     % number of the unique fft points
k = 1:NumUniquePts;         % vector with frequency indexes 
% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of 1/(f^2), 
% i.e. the amplitudes are proportional to 1/f
X = X(1:NumUniquePts);      
X = X./(k.^(beta/2));
% prepare the right half of the spectrum - a conjugate copy of the left
% one except the DC component and the Nyquist component - they are unique,
% and reconstruct the whole spectrum
X = [X conj(X(end-1:-1:2))];
figure(1) 

plot(log(k), log(X(1:NumUniquePts)))
% IFFT
y = real(ifft(X));
% ensure that the length of y is N
y = y(1, 1:N);
% form the noise matrix and ensure unity standard 
% % deviation and zero mean value (columnwise)
% y = reshape(y, [m, n]);
% y = bsxfun(@minus, y, mean(y));
% y = bsxfun(@rdivide, y, std(y));

x=y;
% calculate the noise PSD
winlen = 2*fs;
window = hanning(winlen, 'periodic');
noverlap = winlen/2;
nfft = winlen;
[Pxx, f] = pwelch(x, window, noverlap, nfft, fs, 'onesided');
PxxdB = log10(Pxx);
% plot the noise PSD
figure(2)
semilogx(f, PxxdB, 'r', 'LineWidth', 1.5)
grid on
xlim([1 max(f)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Frequency, Hz')
ylabel('Magnitude, dBV^{2}/Hz')
title('Power Spectral Density of the Noise Signal')