clear all;
close all;

addpath("/home/manny/Documents/Education/11.DSPOctave/functions");
%Exercise 12.2 DSSS Transmission of QPSK
%In this case we apply a Barker code because of its nice spreading property as a short code.

%This program provides simulation for DS-CDMA signaling using coherent QAM detection
%To illustrate the CDMA spreading effect, a single user is spread by the PN sequence of different lengths
%Jamming is added as a narrowband; Changing spreading gain Lc

spreadCode = [ 1 1 1 -1 -1 -1 1 -1 -1 1 -1]';

Ldata = 20000; %The data length must be divisible by 8
Lc = 11; %Spreading factor vs data generate, this can also be changed to 7

%Generate the QPSK modulation symbols
data_sym = 2*round(rand(Ldata, 1))-1 + j*(2*round(rand(Ldata, 1)-1));
jam_data = 2*round(rand(Ldata, 1))-1 + j*(2*round(rand(Ldata, 1)-1));

%Spreading the signals
x_in=kron(data_sym, spreadCode);

%The signal power of the channel input is 2*Lc
%The jamming power is relative
SIR=10;
Pj=2*Lc/(10^(SIR/10));%Power of the jammer

%To generate the AWGN
noiseq=randn(Ldata*Lc,1);+j*randn(Ldata*Lc,1); %The power is 2

%Adding jamming sinusoid sampling frequency, fc = Lc
jam_mod=kron(jam_data,ones(Lc,1)); clear jam_data;
jammer = sqrt(Pj/2)*jam_mod.*exp(j*2*pi*0.12*(1:Ldata*Lc)).'; %fj/fc=0.12, this rotates that amount
clear jam_mod;
% x_in=[1:Ldata];
[P,x]=pwelch(x_in,[],[],[4096],Lc,'twosided');

figure(1);
subplot(2,1,1);semilogy(x-Lc/2, fftshift(P));axis([-Lc/2 Lc/2 1.e-2 1.e2], grid);
xfont=xlabel('frequency (in unit of 1/T_s)');set(xfont, 'Fontsize',11);
yfont=ylabel('CDMA signal PSD');set(yfont, 'Fontsize',11);

[P,x]=pwelch(jammer+x_in,[],[],[4096],Lc,'twosided');
subplot(2,1,2);semilogy(x-Lc/2, fftshift(P));axis([-Lc/2 Lc/2 1.e-2 1.e2], grid);
xfont=xlabel('frequency (in unit of 1/T_s)');set(xfont, 'Fontsize',11);
yfont=ylabel('CDMA signal + Narrowband Jammer PSD');set(yfont, 'Fontsize',11);



BER=[];
BER_az=[];

for i=1:10
  Eb2N(i)=(i-1);
  Eb2N_num=10^(Eb2N(i)/10);
  Var_n=Lc/(2*Eb2N_num);
  signois=sqrt(Var_n);
  awgnois=signois*noiseq;
  %Add noise to signals at the channel output
  y_out=x_in+awgnois+jammer;
  Y_out=reshape(y_out,Lc,Ldata).'; clear y_out awgnois;

  %Despread first
  z_out=Y_out*spreadCode;

  %Decision based on the sign of the samples
  dec1=sign(real(z_out))+j*sign(imag(z_out));
  %Now compare against the original data to compute BER
  BER=[BER; sum([real(data_sym)~=real(dec1);imag(data_sym)~=imag(dec1)])/(2*Ldata)]; %analytical
  BER_az=[BER_az;0.5*erfc(sqrt(Eb2N_num))];
endfor
figure();figber=semilogy(Eb2N,BER_az,'k-',Eb2N,BER,'k-o');
legend('No Jamming','Narrowband Jamming (-10 dB)');set(figber,'LineWidth',2);
xlabel('E_b/N (dB)');ylabel('Bit Error Rate');title('DSSS (CDMA) with spreading gain = 11');
