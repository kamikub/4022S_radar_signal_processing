
clear all;
close all;

fs=20;
f=4;
N=16;
t=(0:1: (N-1))*1/fs;
y = sin(2*pi*f*t);

%% Task 1:
figure();
subplot(2,1,1);
plot(t*1000,y);
xlabel("t (ms)");
ylabel("y");
title("Time Domain:");
subplot(2,1,2);
F=fs*(0:(N-1))/(N-1);
Y=20*log(abs(fft(y)));
plot(F,Y);
xlabel("f (Hz)");
ylabel("Y (dB)");
title("Frequency Domain:");
grid on;                            % Dr Abdul Gaffar

%% Task: 2
figure();
y2=(exp(j*2*pi*f*t));
subplot(2,1,1);
plot(t*1000,abs(y2),t*1000, angle(y2));
xlabel("t (ms)");
ylabel("y");
title("Time Domain:");
subplot(2,1,2);
% F=fs*(0:(N-1))/(N-1);

F=(-N/2:1:(N/2-1))*fs/N;       % Dr Abdul Gaffar

Y=20*log(abs(fftshift(fft(y2))));
plot(F,Y);
xlabel("f (Hz)");
ylabel("Y (dB)");
title("Frequency Domain:");
grid on;                            % Dr Abdul Gaffar


%% Task 3:

figure();
f2=-6;
% y3= cos(2*pi*f*t+0.5)+j*sin(2*pi*f2*t);
y3= exp(j*2*pi*f*t+0.5)+ exp(j*2*pi*f2*t);  % Dr Abdul Gaffar

subplot(2,1,1);
plot(t*1000,abs(y3),t*1000, angle(y3));
xlabel("t (ms)");
ylabel("y");
title("Time Domain:");
subplot(2,1,2);
%F=fs*(0:(N-1))/(N-1);

F=(-N/2:1:(N/2-1))*fs/N;       % Dr Abdul Gaffar

Y=20*log(abs(fftshift(fft(y3))));  % Dr Abdul Gaffar
plot(F,Y);
xlabel("f (Hz)");
ylabel("Y (dB)");
title("Frequency Domain:");
grid on;                            % Dr Abdul Gaffar

%% Task 4
figure();
f2=10;
y4= exp(j*2*pi*f*t+0.5)+ exp(j*2*pi*f2*t);  % Dr Abdul Gaffar
subplot(2,1,1);
plot(t*1000,abs(y4),t*1000, angle(y4));
xlabel("t (ms)");
ylabel("y");
title("Time Domain:");
subplot(2,1,2);
% F=fs*(0:(N-1))/(N-1);
F=(-N/2:1:(N/2-1))*fs/N;           % Dr Abdul Gaffar
Y=20*log(abs(fftshift(fft(y4))));  % Dr Abdul Gaffar
plot(F,Y);
xlabel("f (Hz)");
ylabel("Y (dB)");
title("Frequency Domain:");
grid on;                            % Dr Abdul Gaffar
