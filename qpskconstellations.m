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
    symbolToAdd += frame(k+(q-1))*2^(bitsPerSymbol-q);
  endfor
  symbolToAdd+=1; %Because I am using based-1
  dataSymbols = [dataSymbols symbolToAdd*ones(1,number) ];
endfor

dataEncoded = symbolMap(dataSymbols);

% *********************************************************************
% Now lets modulate the signal

T = 1;
Fs = 200; % Sampling frequency
Ts = 1/Fs; %Sampling Period
t=0.005:Ts:length(dataEncoded)*Ts;
f=5;%Frequency of the carrier
%Here we generate the modulated signal by multiplying it with
%carrier (basis function)
% Modulated=abs(dataEncoded).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
% Modulated=exp(j*pi*arg(dataEncoded)).*(sqrt(2/T)*cos(2*pi*f*t));
Modulated=(sqrt(2/T)*cos(2*pi*f*t + arg(dataEncoded)));%.*exp(j*pi*arg(dataEncoded));

L = length(dataEncoded);
freq = Fs*(0:(L/2))/L;
sig_fft = fft(Modulated);
P2 = abs(sig_fft/L); %Computes the two-sided spectrum
P1 = P2(1:L/2+1); %Computes the single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1); %Computes the even-valued signal length L

subplot(4,2,1:2);
maxValue = max(Modulated);
plot(Modulated,".-");axis([0 number*4 -maxValue maxValue], grid);xlabel('Time (seconds)');ylabel('Amplitude (volts)');title('QPSK Modulated signal');

strFrame = ["Frame UW: " mat2str(frame(1:length(uw))) ];
subplot(4,2,3:4);
plotFrame = kron(frame, ones(1,(number/bitsPerSymbol)));
stem(plotFrame,".-");axis([-0.1 number*4 -0.1 1.1], grid);xlabel('Bit Position');ylabel('Value');title(strFrame);
% polar(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');

subplot(4,2,5);
polar(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');
subplot(4,2,6);
maxValue = max(P1);
plot(freq, P1);axis([-f*1 f*3 0 maxValue], grid);xlabel("f\(Hz\)");ylabel("|P1(f)|");
xlabel('frequency (Hz)');ylabel('Amplitude');title('Power Spectrum');
% plot(dataEncoded,".");axis(grid);xlabel('Re');ylabel('Im');title('QPSK Modulated signal');

transmittedSymbols = dataSymbols(1:number:number*5);
transmittedSamples = dataEncoded(1:number:end);
transmittedPhases = arg(transmittedSamples)*180/pi;

printf("The first transmitted phases are: %s = %d; %s = %d; %s = %d; %s = %d; %s = %d\n",
    dec2bin(transmittedSymbols(1)-1), transmittedPhases(1), dec2bin(transmittedSymbols(2)-1), transmittedPhases(2),
    dec2bin(transmittedSymbols(3)-1), transmittedPhases(3), dec2bin(transmittedSymbols(4)-1), transmittedPhases(4),
    dec2bin(transmittedSymbols(5)-1), transmittedPhases(5));

countPhase00 = length(find(transmittedPhases == 45));
countPhase01 = length(find(transmittedPhases == 135));
transmittedPhasesRemovePositive = transmittedPhases(find(transmittedPhases != 135));
transmittedPhasesRemovePositive = transmittedPhasesRemovePositive(find(transmittedPhasesRemovePositive != 45));
transmittedPhasesRemovePositive(1:end) *= -1.0;
countPhase10 = length(find(transmittedPhasesRemovePositive == 135));
countPhase11 = length(find(transmittedPhasesRemovePositive == 45));

printf("The total count of transmitted phases are: 00=>(45) = %d; 01=>(135) = %d; 10=>(-135) = %d; 11=>(-45) = %d\n",
    countPhase00, countPhase01, countPhase10, countPhase11);


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
