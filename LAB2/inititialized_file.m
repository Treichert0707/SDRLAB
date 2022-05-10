function SimParams = inititialized_file
%   Copyright 2017 The MathWorks, Inc.
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
SimParams.Rsym = 0.4e6;             % Symbol rate in Hertz
SimParams.ModulationOrder = 4;      % QPSK alphabet size
SimParams.Interpolation = 2;        % Interpolation factor
SimParams.Decimation = 1;           % Decimation factor
SimParams.Tsym = 1/SimParams.Rsym;  % Symbol time in sec
SimParams.Fs   = SimParams.Rsym * SimParams.Interpolation; % Sample rate

%% Frame Specifications
% [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ... | 'Hello world 099\n'];
SimParams.BarkerCode      = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];     % Bipolar Barker Code
SimParams.BarkerLength    = length(SimParams.BarkerCode);
SimParams.HeaderLength    = SimParams.BarkerLength * 2;                   % Duplicate 2 Barker codes to be as a header
SimParams.Message         = message_bits;
SimParams.MessageLength   = length(message_bits);                % 'Hello world 000\n'...
SimParams.NumberOfMessage = 1;                                          % Number of messages in a frame
SimParams.PayloadLength   = length(message_bits); % 7 bits per characters
SimParams.FrameSize       = (SimParams.HeaderLength + SimParams.PayloadLength) ...
    / log2(SimParams.ModulationOrder);                                    % Frame size in symbols
SimParams.FrameTime       = SimParams.Tsym*SimParams.FrameSize;

%% Rx parameters
SimParams.RolloffFactor     = 0.5;                      % Rolloff Factor of Raised Cosine Filter
SimParams.ScramblerBase     = 2;
SimParams.ScramblerPolynomial           = [1 1 1 0 1];
SimParams.ScramblerInitialConditions    = [0 0 0 0];
SimParams.RaisedCosineFilterSpan = 10;                  % Filter span of Raised Cosine Tx Rx filters (in symbols)
SimParams.DesiredPower                  = 2;            % AGC desired output power (in watts)
SimParams.AveragingLength               = 50;           % AGC averaging length
SimParams.MaxPowerGain                  = 60;           % AGC maximum output power gain
SimParams.MaximumFrequencyOffset        = 6e3;
% Look into model for details for details of PLL parameter choice. 
% Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice.
K = 1;
A = 1/sqrt(2);
SimParams.PhaseRecoveryLoopBandwidth    = 0.01;         % Normalized loop bandwidth for fine frequency compensation
SimParams.PhaseRecoveryDampingFactor    = 1;            % Damping Factor for fine frequency compensation
SimParams.TimingRecoveryLoopBandwidth   = 0.01;         % Normalized loop bandwidth for timing recovery
SimParams.TimingRecoveryDampingFactor   = 1;            % Damping Factor for timing recovery
% K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM),
% QPSK could be treated as two individual binary PAM,
% 2.7 is for raised cosine filter with roll-off factor 0.5
SimParams.TimingErrorDetectorGain       = 2.7*2*K*A^2+2.7*2*K*A^2;
SimParams.PreambleDetectorThreshold     = 0.8;

%% Message generation and BER calculation parameters
%msgSet = zeros(100 * SimParams.MessageLength, 1);
%%for msgCnt = 0 : 99
%    msgSet(msgCnt * SimParams.MessageLength + (1 : SimParams.MessageLength)) = ...
 %       sprintf('%s %03d\n', SimParams.Message, msgCnt);
%end
%bits = de2bi(msgSet, 7, 'left-msb')';
SimParams.MessageBits = message_bits;

% For BER calculation masks
SimParams.BerMask = zeros(length(SimParams.MessageBits), 1);
% Pluto receiver parameters
SimParams.PlutoCenterFrequency = 2.2e9 - 1.7e4;
SimParams.PlutoGain                 = 15;
SimParams.PlutoFrontEndSampleRate   = SimParams.Fs;
SimParams.PlutoFrameLength          = SimParams.Interpolation * SimParams.FrameSize;

% Experiment parameters
SimParams.PlutoFrameTime = SimParams.PlutoFrameLength / SimParams.PlutoFrontEndSampleRate;
SimParams.StopTime = 4.8;
