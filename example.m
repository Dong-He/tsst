close all;clear;clc;
% This code is to test the tsst
% Creat the test signal
N=1024;
fs=1024;
t=(0:N-1)/fs;
sigma=0.00001;
t0=0.5;
x = (pi*sigma^2)^(-0.25).* exp(-((t-t0)/sigma).^2/2);
x = x(:);
% time waveform
figure
plot(t,x);xlabel('Time/s');title('Time Waveform')
% analysis parameters
SIGMA=0.01;
% perform tsst
[STFT,Vx,time,freqr]=time_synsq_stft_fw(x,fs,SIGMA,1);
% figure output
figure
subplot(121)
imagesc(time,freqr,abs(STFT));axis xy
xlabel('Time/s');ylabel('Frequency/Hz');title('STFT')
subplot(122)
imagesc(time,freqr,abs(Vx));axis xy
xlabel('Time/s');ylabel('Frequency/Hz');title('TSST')
