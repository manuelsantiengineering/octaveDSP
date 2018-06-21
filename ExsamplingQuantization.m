clear all;
close all;

Fd = 500; %Sampling rate of 500Hz
td = 1/Fd;
t = [0:td:1];
f1 = 1;
f2 = 3;
xsig = sin(2*pi*f1*t)-sin(2*pi*f2*t);
Lsig = length(xsig);


sr = 50; %Sampling rate of 50Hz
ts = 1/sr;
Nfactor = sr/td;
%Using a 16 level uniform quantizer
[ s_out, sq_out, sqh_out, Delta, SQNR ] = functsampandquant(xsig,16, td, ts );
Lenfft = 2^ceil(log2(Lsig)+1); %TO make it a power of 2
Fmax=0.5*Fd;
Faxis=linspace(-Fmax,Fmax,Lenfft);
Xsig=fftshift(fft(xsig,Lenfft));
S_out=fftshift(fft(s_out,Lenfft));

figure(1);
subplot(3,1,1);sfig1a=plot(t,xsig,"k");axis(grid);
hold on; sfig1b=plot(t,s_out(1:Lsig),"b");hold off;legend("Original","Sampled");
set(sfig1a,"linewidth",2);set(sfig1b,"linewidth",2.); xlabel("time(sec)");
title("Signal and its uniform samples");
subplot(3,1,2);sfig1c=plot(Faxis,abs(Xsig),"k");axis(grid);
set(sfig1c,"linewidth",1);xlabel("freq(Hz)");
subplot(3,1,3);sfig1d=plot(Faxis,abs(S_out),"b");axis(grid);
set(sfig1c,"linewidth",1);xlabel("freq(Hz)");
