clc;
clear all;
close all;
[y_n, Fs] = audioread('sf3_n0H.wav');
[y_clean, Fs] = audioread('sf3_cln.wav');

dt = 1/Fs;
t = 0:dt:(length(y_n)*dt)-dt;
%% Generate figures of the noisey and clean signal
L = length(y_n);

Y_n = fft(y_n);
P2 = abs(Y_n/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Y_clean = fft(y_clean);
O2 = abs(Y_clean/L);
O1 = O2(1:L/2+1);
O1(2:end-1) = 2*O1(2:end-1);

f1 = Fs*(0:(L/2))/L;
f2 = Fs*(0:(L/2))/L;
subplot(2,2,1)
plot(f1,P1) 
title('Single-Sided Amplitude Spectrum of Distorted Audio')
xlabel('f (Hz)')
ylabel('|Ydistorted(f)|')
axis([0 8000 0 max(P1)])
hold on
subplot(2,2,2)
spectrogram(y_n)
subplot(2,2,3)
plot(f2,O1) 
title('Single-Sided Amplitude Spectrum of Clean Audio')
xlabel('f (Hz)')
ylabel('|Yclean(f)|')
axis([0 8000 0 max(P1)])
subplot(2,2,4)
spectrogram(y_clean)

%% Configure the adaptive filter
filter_length = 16; % filter order
x_n = y_n; % input signal
signal_length = length(x_n);
y = zeros(length(y_n),1);
e = zeros(length(y_n),1);
w = zeros(1,filter_length);
x = zeros(1,filter_length);
mu = zeros(1,signal_length);

for i = 1:signal_length
    x(1) = x_n(i);
    y(i) = w*x';
    e(i) = y_clean(i) - y(i);
    mu(i) = 2/max(eig((x'*x)));
    w = w + mu(i)*e(i)*x;   
    x(2:end) = x(1:end-1);
end


figure
subplot(4,1,1)
plot(t,y_n); xlabel('Seconds'); ylabel('Amplitude');
hold on
subplot(4,1,2)
plot(psd(spectrum.periodogram,y_n,'Fs',Fs,'NFFT',length(y_n)));
subplot(4,1,3)
plot(t,y); xlabel('Seconds'); ylabel('Amplitude');
subplot(4,1,4)
plot(psd(spectrum.periodogram,y,'Fs',Fs,'NFFT',length(y)));


figure;
subplot(3,1,1)
plot(t,y_n); xlabel('Seconds'); ylabel('Amplitude');
title('Noisey Signal')
hold on;
subplot(3,1,2)
plot(t,y); xlabel('Seconds'); ylabel('Amplitude');
title('Adapted Signal')
subplot(3,1,3)
plot(t,y_clean); xlabel('Seconds'); ylabel('Amplitude');
title('Clean Signal')

figure;
zplane(w,1)

figure;
plot(200:signal_length,mu(200:end))
xlabel('Signal Length')
ylabel('Value of Mu')

fvtool(w)
phasez(w,1,2048)


