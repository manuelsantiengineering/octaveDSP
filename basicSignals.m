clear all;
% close all;
%

% sr = 20;
% del = 1/sr;
% f = 1;
% T0 = 1/f;
% sr = 1;
% n = 0:100;
% A = 1;
% x = A*sin(2*pi*n*del/T0);

Fs = 100000; %Sampling Freq
Ts = 1/Fs; %Sampling Period
L = 25000; %Length of signal
t = (0:L-1)*Ts; %Time vector
A = 1;

f1 = 10;
T = (1/f1);
%
% sig = A*sin(2*pi*f1*t);
%
% freq = Fs*(0:(L/2))/L;
% sig_fft = fft(sig);
% P2 = abs(sig_fft/L); %Computes the two-sided spectrum
% P1 = P2(1:L/2+1); %Computes the single-sided spectrum
% P1(2:end-1) = 2*P1(2:end-1); %Computes the even-valued signal length L
%
% subplot(2,2,1);
% % stem(sig);axis([0 T*2 -A A], grid);xlabel("t\(sec\)");ylabel("sin(2*pi*f*t)");
% stem(t, sig);axis([0 2*T -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,2);
% plot(t, sig);axis([0 2*T -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,3);
% plot(freq, P1);axis([-f1*2 f1*2 0 1.1], grid);xlabel("f\(Hz\)");ylabel("|P1(f)|");
% subplot(2,2,4);
% plot3(sig);axis([0 1000*T -A A -A A], grid);xlabel("t\(ms\)");ylabel("A");
% % f = 2; %Freq tells how may cycles per second
% % cycle = 2*pi; %For or a sinusoid, a full cycle occurs at 2*pi


sig = 0;
for freq = 1:2:200
  sig += sin(2*pi*freq*t);
endfor

subplot(2,2,1);
% stem(sig);axis([0 T*2 -A A], grid);xlabel("t\(sec\)");ylabel("sin(2*pi*f*t)");
stem(t, sig);axis([0 1 -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
subplot(2,2,2);
plot(t, sig);axis([0 2*T -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
subplot(2,2,3);
plot(freq, P1);axis([-f1*2 f1*2 0 1.1], grid);xlabel("f\(Hz\)");ylabel("|P1(f)|");
subplot(2,2,4);
plot3(sig);axis([0 1000*T -A A -A A], grid);xlabel("t\(ms\)");ylabel("A");

% T = 1/f;
% sr = 5;
% ts = 1/sr;
%
% samples = [0:ts:cycle*T];
% t = samples/cycle;
% A = 1;
%
% sig = A*exp(pi*j*[0:0.1:2*pi]);
% rot_angle = 180;
% rot_rad = rot_angle*(pi/180);
% % sig_rot = sig.*exp(pi*j*rot_rad*[0:length(sig)-1]);
% sig_rot = sig.*exp(pi*j*rot_rad);
% % sig_cos = cos(samples*f);
% % sig_sin = j*sin(samples*f);
% % sig_sin_rot = sig_sin.*exp(2*pi*j*[0:length(sig_sin)-1]);
% % sig = sig_cos + j*sig_sin;
% % sig = exp(2*pi*j*[0:0.2:2*pi]);
% % sig = A*sin(samples*f);
% % subplot(1,2,1);
% % scatter(sig_cos, sig_sin);axis([-1.1 1.1 -1.1 1.1]*A, grid);
% subplot(1,2,1);
% % plot(t, sig, "-x");axis(grid);%xlabel("cos(");ylabel("sin(x)");
% plot(sig, "rx");axis(grid);%xlabel("cos(");ylabel("sin(x)");
% % plot(sig, "-x");
% subplot(1,2,2);
% plot(sig_rot, "bx");axis(grid);%xlabel("cos(");ylabel("sin(x)");
