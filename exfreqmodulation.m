clear all;
close all;

addpath("/home/manny/Documents/Education/11.DSPOctave/functions");

ts = 1.e-4;

t = -0.04:ts:0.04;
Ta=0.01;
m_sig=triangl((t+0.01)/Ta)-triangl((t-0.01)/Ta);
Lfft = length(t); Lfft=2^ceil(log2(Lfft));

M_fre=fftshift(fft(m_sig, Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=100; %Bandwidth of the signal (Hz)

%Display a simple lowpass filter width bandwidth B_m Hz
h=fir1(80, [B_m*ts]);

kf=160*pi;
m_intg=kf*ts*cumsum(m_sig);
f=300;
s_fm=cos(2*pi*f*t+m_intg);
s_pm=cos(2*pi*f*t+pi*m_sig);
S_fm=fftshift(fft(s_fm,Lfft));
S_pm=fftshift(fft(s_pm,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_fmdem=diff([s_fm(1) s_fm])/ts/kf;
s_fmrec=s_fmdem.*(s_fmdem>0);
s_dec=filter(h,1,s_fmrec);

%Demodulation
%Using an ideal LPF with bandwidth 200 Hz
Trange1=[-0.04 0.04 -1.2 1.2];

figure(1);
subplot(2,2,1);m1=plot(t,m_sig);axis(grid,Trange1);set(m1,"Linewidth",2);
xlabel('{\it t} (sec)');ylabel('{\it m} ({\it t})');title('Message Signal');

subplot(2,2,2);m2=plot(t,s_dec);axis(grid);set(m2,"Linewidth",2);
xlabel('{\it t} (sec)');ylabel('{\it m}_d ({\it t})');title('Demodulated FM Signal');

subplot(2,2,3);td1=plot(t,s_fm);axis(grid,Trange1);set(td1,"Linewidth",2);
xlabel('{\it t} (sec)');ylabel('{\it s}_{\rm FM} ({\it t})');title('FM Signal');

subplot(2,2,4);td2=plot(t,s_pm);axis(grid,Trange1);set(td2,"Linewidth",2);
xlabel('{\it t} (sec)');ylabel('{\it s}_{\rm PM} ({\it t})');title('PM Signal');

figure(2);
subplot(2,2,1);fp1=plot(t,s_fmdem);axis(grid);set(fp1,"Linewidth",2);
xlabel('{\it t} (sec)');ylabel('{\it s}_{\rm FM} ({\it t})');title('FM Derivative');

subplot(2,2,2);fp2=plot(t,s_pm);axis(grid,Trange1);set(fp2,"Linewidth",2);
xlabel('{\it t} (sec)');title('Rectified FM Derivative');

subplot(2,2,3);fp1=plot(freqs,abs(S_fm));axis(grid);set(fp1,"Linewidth",2);
xlabel('{\it f} (Hz)');ylabel('{\it S}_{\rm FM} ({\it f})');title('FM Amplitude Spectrum');

subplot(2,2,4);fp2=plot(freqs,abs(S_pm));axis(grid);set(fp2,"Linewidth",2);
xlabel('{\it f} (Hz)');ylabel('{\it S}_{\rm PM} ({\it f})');title('PM Amplitude Spectrum');
