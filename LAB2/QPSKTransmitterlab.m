Fs = 8000;
M = 4;
samplingFactor = 1;
Mu = 2.^8-1;
samples = [1,Fs];
[y, Fs] = audioread('audio_file.wav',samples);
soundData = resample(y,M,1); % This upsamples the signal by M
downSampledSignal = downsample(soundData,samplingFactor); % This downsamples the signal by samplingfactor
maximumValueofSignal = max(downSampledSignal);
compandedSignal = compand(downSampledSignal,Mu,maximumValueofSignal,'mu/compressor');
maxofCompandedSignal = max(soundData);
minofCompandedSignal = min(soundData);
stepSize = (maxofCompandedSignal - minofCompandedSignal)/Mu;
partition =[minofCompandedSignal:stepSize:maxofCompandedSignal];
codeBook = [(minofCompandedSignal-stepSize):stepSize:maxofCompandedSignal];
[index,quants] = quantiz(compandedSignal,partition,codeBook);
sourceCode = [];
sourceCode = [sourceCode,dec2bin(index,8)];
souceCodeTranspose = sourceCode';
sourceCodeStream = reshape(souceCodeTranspose,1,[]);
symbolValue = []; % initilization
counter = 0; %initilization

for bit = 1:length(sourceCodeStream)/2
symbolValue = [symbolValue;sourceCodeStream(bit+counter:bit+counter+1)];
counter = counter +1;
end
decimalSymbolValue = bin2dec(symbolValue);
modulatedSignal = pskmod(decimalSymbolValue,M);
for bit_maker = 1:256000
bit_value_double(bit_maker) = str2double(sourceCodeStream(bit_maker));
end
message_bits = bit_value_double';

%% General simulation parameters
SimParams.Rsym = 0.4e6; % Symbol rate in Hertz
SimParams.ModulationOrder = 4; % QPSK alphabet size
SimParams.Interpolation = 2; % Interpolation factor
SimParams.Decimation = 1; % Decimation factor
SimParams.Tsym = 1/SimParams.Rsym; % Symbol time in sec
SimParams.Fs = SimParams.Rsym * SimParams.Interpolation; % Sample rate
%% Frame Specifications (correct these as you need to)
% [BarkerCode*2 |Message_Bits];
SimParams.BarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; % BipolarBarker Code
SimParams.BarkerLength = length(SimParams.BarkerCode);
SimParams.HeaderLength = SimParams.BarkerLength * 2; % Duplicate 2 Barker codes to be as a header
%SimParams.Message = 'HELLO WORLD'; % No longer used
SimParams.MessageLength = length(message_bits); % the length of the new message bits
SimParams.NumberOfMessage = 1;
SimParams.PayloadLength = length(message_bits);
SimParams.FrameSize = (SimParams.HeaderLength + SimParams.PayloadLength) ...
 / log2(SimParams.ModulationOrder); % Framesize in symbols
SimParams.FrameTime = SimParams.Tsym*SimParams.FrameSize;
SimParams.MessageBits = message_bits;
%% RRC Parameters
SimParams.RolloffFactor = 0.5; % Rolloff Factor of Raised Cosine Filter
SimParams.ScramblerBase = 2;
SimParams.ScramblerPolynomial = [1 1 1 0 1];
SimParams.ScramblerInitialConditions = [0 0 0 0];
SimParams.RaisedCosineFilterSpan = 10; % Filter span of Raised Cosine Tx Rx filters (in symbols)
% Pluto transmitter parameters
SimParams.PlutoCenterFrequency = 2.2e9;
SimParams.PlutoGain = 0;
SimParams.PlutoFrontEndSampleRate = SimParams.Fs;
SimParams.PlutoFrameLength = SimParams.Interpolation * SimParams.FrameSize;
% Simulation Parameters
SimParams.FrameTime = SimParams.PlutoFrameLength/SimParams.PlutoFrontEndSampleRate;
SimParams.StopTime = 1000;
SimParams.Address = 'usb:0';

%%
runPlutoradioQPSKTransmitterOK(SimParams);





