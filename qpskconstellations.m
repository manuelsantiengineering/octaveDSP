clear all;
% close all;

% *********************************************************************
% This part is just mappingthe constellation and assigning the numbering

bitsPerSymbol = 2; %Meaning it is QPSK

map = exp(pi*j*[1:2:7]/4);
numbering = [1 2 3 4]; %Octave is base 1

% map16QAM = [];
x = -3 + j*[-3:2:3];
map16QAM = (x')+[0:2:6];


% map = exp(pi*j*[0:1]);
% numbering = [1 2]; %Octave is base 1

symbolMap = map(numbering);

%This part is just used to understand how symbols are mapped.
angleSymbolMap = arg(symbolMap)*180/pi;
symbols = [0b00 0b01 0b10 0b11]; % use dec2bin()
symbolsArranged = symbols(numbering)+1;
symbols = dec2bin(symbolsArranged-1)';
printf("The symbol mapping is:\t %s = %d; %s = %d; %s = %d; %s = %d\n",
    symbols(1:2), angleSymbolMap(1), symbols(3:4), angleSymbolMap(2),
    symbols(5:6), angleSymbolMap(3), symbols(7:8), angleSymbolMap(4));
% *********************************************************************
% This part is just constructing the bit frames
% Frame = [10 bits UW][8 bits Length][N bits Data]
%The length bits are interpreted in their decimal value and using 8 bits, so the maximum is 255
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
f=5;%Frequency of the carrier
%Here we generate the modulated signal by multiplying it with
%carrier (basis function)
% Modulated=abs(dataEncoded).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
Modulated=(sqrt(2/T)*cos(2*pi*f*t + arg(dataEncoded)));%.*exp(j*pi*arg(dataEncoded));
maxValue = max(Modulated);
subplot(4,2,1:2);
plot(Modulated,".-");axis([0 number*4 -maxValue maxValue], grid);xlabel('Time (seconds)');ylabel('Amplitude (volts)');title('QPSK Modulated signal');

subplot(4,2,3:4);
plotFrame = kron(frame, ones(1,(number/bitsPerSymbol)));
stem(plotFrame,".-");axis([-0.1 number*4 -0.1 1.1], grid);xlabel('Bit Position');ylabel('Value');title('Frame');
% polar(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');

subplot(4,2,5);
polar(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');
subplot(4,2,6);
plot(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');

receivedSymbols = dataSymbols(1:number:number*5);
receivedSamples = dataEncoded(1:number:number*5);
receivedPhases = arg(receivedSamples)*180/pi;

printf("The received phases are: %s = %d; %s = %d; %s = %d; %s = %d; %s = %d\n",
    dec2bin(receivedSymbols(1)-1), receivedPhases(1), dec2bin(receivedSymbols(2)-1), receivedPhases(2),
    dec2bin(receivedSymbols(3)-1), receivedPhases(3), dec2bin(receivedSymbols(4)-1), receivedPhases(4),
    dec2bin(receivedSymbols(5)-1), receivedPhases(5));

% Fs = 200; %Sampling Freq
% Ts = 1/Fs; %Sampling Period
% L = length(dataEncoded); %Length of signal
% t1 = (0.005:L-1)*Ts; %Time vector
% A = 1;
%
% f1 = 10;
% T = (1/f1);

% *********************************************************************
% So the signal is already modulated, and "transmitted"
% It is time to demodulate it



%We begin demodulation by multiplying the received signal again with
%the carrier (basis function)
demodulated=Modulated.*(sqrt(2/T)*cos(2*pi*f*t));
%Here we perform the integration over time period T using trapz
%Integrator is an important part of correlator receiver used here
y=[];
for i=1:number:size(demodulated)(2)
 y=[y trapz(t(i:i+(number-1)),demodulated(i:i+(number-1)))];
end
received=y>0;

maxValue = max(demodulated);
subplot(4,2,7:8);
plot(demodulated,".-");axis([0 number*4 -maxValue maxValue], grid);xlabel('Time (seconds)');ylabel('Amplitude (volts)');title('QPSK Demodulated signal');
