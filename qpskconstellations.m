clear all;
% close all;

% *********************************************************************
% This part is just mappingthe constellation and assigning the numbering

bitsPerSymbol = 2; %Meaning it is QPSK

map = exp(pi*j*[1:2:7]/4);
numbering = [1 2 3 4]; %Octave is base 1

map16QAM = [];
% map = exp(pi*j*[0:1]);
% numbering = [1 2]; %Octave is base 1

symbolMap = map(numbering)

%This part is just used to understand how symbols are mapped.
angleSymbolMap = arg(symbolMap)*180/pi
symbols = [0b00 0b01 0b10 0b11]; % use dec2bin()
symbolsArranged = symbols(numbering)+1;
symbols = dec2bin(symbolsArranged-1)

% *********************************************************************
% This part is just constructing the bit frames
% Frame = [10 bits UW][8 bits Length][N bits Data]
%The length bits are interpreted in their decimal value, so the maximum is 255
stringToModulate = "Hola";
uw = [1 1 0 0 0 1 1 0 1 1]
dataBits = dec2bin(toascii(stringToModulate));
dataColNum =  size(dataBits)(2); %Used to verify if it is using 7 bits or 8 bits for the hex
dataRowNum = size(dataBits)(1);

dataBits = reshape(dataBits',1,[]);
dataBitsStr = int2str(dataBits);
dataBits = arrayfun( @(dataBits) str2double(dataBitsStr(dataBits)), 1:numel(dataBitsStr) );

bitsToInterleave = [];

if(dataColNum < 8)
  numOfBitsToInterleave = 8-dataColNum;
  for k = 1:numOfBitsToInterleave
    bitsToInterleave = [bitsToInterleave 0];
  endfor

  % dataBits = [bitsToInterleave dataBits];
  % interleavePosition = dataColNum+numOfBitsToInterleave+1;
  dataBitsInterleaved = [];
  for k = 0:(dataRowNum-1)
    startBit = 1+dataColNum*i;
    endBit = startBit + (dataColNum-1);
    dataBitsInterleaved = [ dataBitsInterleaved bitsToInterleave dataBits(startBit:endBit)];
  endfor
  dataBits = dataBitsInterleaved;
endif

dataBitsLength = dec2bin(length(dataBits));
dataBitsStr = int2str(dataBitsLength);
dataBitsLength = arrayfun( @(dataBitsLength) str2double(dataBitsStr(dataBitsLength)), 1:numel(dataBitsStr) );
numOfBitsToInterleave = 8 - length(dataBitsLength); %Used to verify if it is using 7 bits or 8 bits for the hex

bitsToInterleave = [];
for k = 1:numOfBitsToInterleave
  bitsToInterleave = [bitsToInterleave 0];
endfor
dataBitsLength = [bitsToInterleave dataBitsLength];

frame = [uw dataBitsLength dataBits];
% *********************************************************************
% This part is where the bits are converted to analog values

number = 200;
dataSymbols = [];
for k = 1:bitsPerSymbol:length(frame)
  symbolToAdd = 0;
  for q = 1:(bitsPerSymbol)
    symbolToAdd += frame(q)*2^(bitsPerSymbol-q);
  endfor
  symbolToAdd+=1; %Because I am using based-1
  dataSymbols = [dataSymbols (((frame(k)*2^1)+frame(k+1))+1)*ones(1,number) ];
  % dataSymbols = [dataSymbols symbolToAdd*ones(1,number) ];
endfor

dataEncoded = symbolMap(dataSymbols);

% *********************************************************************
% Now lets modulate the signal

T = 1;
timeStep = 0.005;
t=0.005:timeStep:length(dataEncoded)*timeStep;
%Frequency of the carrier
f=5;
%Here we generate the modulated signal by multiplying it with
%carrier (basis function)
% Modulated=abs(dataEncoded).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
Modulated=(sqrt(2/T)*cos(2*pi*f*t + arg(dataEncoded)));%.*exp(j*pi*arg(dataEncoded));
maxValue = max(Modulated);
subplot(2,2,1:2);
plot(Modulated,".-");axis([0 number*4 -maxValue maxValue], grid);
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('QPSK Modulated signal');

subplot(2,2,3);
polar(dataEncoded,".");axis(grid);
xlabel('Re');
ylabel('Im');
title('QPSK Modulated signal');

subplot(2,2,4);
plot(dataEncoded,".");axis(grid);
xlabel('Re');
ylabel('Im');
title('QPSK Modulated signal');

dataSymbols(1:number:number*8)
dataEncoded(1:number:number*8)
arg(dataEncoded(1:number:number*8))*180/pi



% Fs = 10000; %Sampling Freq
% Ts = 1/Fs; %Sampling Period
% L = 25000; %Length of signal
% t = (0:L-1)*Ts; %Time vector
% A = 1;
%
% f1 = 10;
% T = (1/f1);





% polar(dataAnalog, "r.-");title("Original Analog Signal");
% plot(dataAnalog, "r.-");axis(grid);title("Original Analog Signal");

% subplot(2,2,1);
% stem(t, sig);axis([0 1 -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,2);
% plot(t, sig);axis([0 2*T -A A ], grid);xlabel("t\(s\)");ylabel("sin(2*pi*f*t)");
% subplot(2,2,3);
% plot(freq, P1);axis([-f1*2 f1*2 0 1.1], grid);xlabel("f\(Hz\)");ylabel("|P1(f)|");
% subplot(2,2,4);
% plot3(sig);axis([0 1000*T -A A -A A], grid);xlabel("t\(ms\)");ylabel("A");
