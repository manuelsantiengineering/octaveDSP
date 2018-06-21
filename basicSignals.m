clear all;
close all;

%
% Fs = 100000; %Sampling Freq
% Ts = 1/Fs; %Sampling Period
% L = 25000; %Length of signal
% t = (0:L-1)*Ts; %Time vector
% A = 1;
%
% f1 = 10;
% T = (1/f1);
% %
% sig = A*sin(2*pi*f1*t);
%
% freq = Fs*(0:(L/2))/L;
% sig_fft = fft(sig);
% P2 = abs(sig_fft/L); %Computes the two-sided spectrum
% P1 = P2(1:L/2+1); %Computes the single-sided spectrum
% P1(2:end-1) = 2*P1(2:end-1); %Computes the even-valued signal length L
%
%
% subplot(2,2,1);
% % stem(sig);axis([0 T*2 -A A], grid);xlabel("t\(sec\)");ylabel("sin(2*pi*f*t)");
% stem(t, sig);axis([0 1 -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,2);
% plot(t, sig);axis([0 2*T -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,3);
% plot(freq, P1);axis([-f1*2 f1*2 0 1.1], grid);xlabel("f\(Hz\)");ylabel("|P1(f)|");
% subplot(2,2,4);
% plot3(sig);axis([0 1000*T -A A -A A], grid);xlabel("t\(ms\)");ylabel("A");

f = 10; %frequency
