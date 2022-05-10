%% Concatenate the Frames into a single column vector and pass these into a buffer_signal upon which we will operate.
%your code for concatenation
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


ConcatenatedSignal = reshape(recorded_signal,[],1);
bufferSignal = ConcatenatedSignal;
obj.pAGC = comm.AGC( ...
 'DesiredOutputPower', prmQPSKReceiver.DesiredPower, ...
 'AveragingLength', prmQPSKReceiver.AveragingLength, ...
 'MaxPowerGain', prmQPSKReceiver.MaxPowerGain);
AGCSignal = obj.pAGC(bufferSignal);
obj.pRxFilter = comm.RaisedCosineReceiveFilter( ...
 'RolloffFactor', prmQPSKReceiver.RolloffFactor, ...
 'FilterSpanInSymbols', prmQPSKReceiver.RaisedCosineFilterSpan, ...
 'InputSamplesPerSymbol', prmQPSKReceiver.Interpolation, ...
 'DecimationFactor', prmQPSKReceiver.Decimation);
RCRxSignal = obj.pRxFilter(AGCSignal);
obj.pCoarseFreqEstimator = comm.CoarseFrequencyCompensator( ...
 'Modulation', 'QPSK', ...
 'Algorithm', 'Correlation-based', ...
 'MaximumFrequencyOffset', prmQPSKReceiver.MaximumFrequencyOffset, ....
 'SampleRate', prmQPSKReceiver.Fs);

 obj.pCoarseFreqCompensator = comm.PhaseFrequencyOffset( ...
 'PhaseOffset', 0, ...
 'FrequencyOffsetSource', 'Input port', ...
 'SampleRate', prmQPSKReceiver.Fs);
obj.pMeanFreqOff = 0;

obj.pCnt = 0;

 obj.pFineFreqCompensator = comm.CarrierSynchronizer( ...
 'Modulation', 'QPSK', ...
 'ModulationPhaseOffset', 'Auto', ...
 'SamplesPerSymbol', prmQPSKReceiver.Interpolation, ...
 'DampingFactor', prmQPSKReceiver.Decimation, ...
 'NormalizedLoopBandwidth', prmQPSKReceiver.PhaseRecoveryLoopBandwidth);
[~, freqOffsetEst] = obj.pCoarseFreqEstimator(RCRxSignal); % Coarse frequency offset estimation
 % average coarse frequency offset estimate, so that carrier
 % sync is able to lock/converge
 freqOffsetEst = (freqOffsetEst + obj.pCnt * obj.pMeanFreqOff)/(obj.pCnt+1);
 obj.pCnt = obj.pCnt + 1; % update state
 obj.pMeanFreqOff = freqOffsetEst;
 coarseCompSignal = obj.pCoarseFreqCompensator(RCRxSignal,...
 -freqOffsetEst); % Coarse frequency compensation
 obj.pTimingRec = comm.SymbolSynchronizer( ...
 'TimingErrorDetector', 'Gardner (non-data-aided)', ...
 'SamplesPerSymbol', prmQPSKReceiver.Interpolation, ...
 'DampingFactor', prmQPSKReceiver.Decimation, ...
 'NormalizedLoopBandwidth',prmQPSKReceiver.PhaseRecoveryLoopBandwidth, ...
 'DetectorGain', 5.4);
 timingRecSignal = obj.pTimingRec(coarseCompSignal); % Symbol timing recovery
 fineCompSignal = obj.pFineFreqCompensator(timingRecSignal); % Fine frequency compensation

 %%
real_val = real(fineCompSignal);
imag_val = imag(fineCompSignal);
scatter(real_val, imag_val)
abs(max(fineCompSignal))
obj.pModulatedHeader = sqrt(2)/2 * (-1-1i) * [+1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1];
obj.pPrbDet = comm.PreambleDetector(obj.pModulatedHeader, ...
 'Input', 'Symbol', ...
 'Threshold', 25.2, ...
 'Detections', 'All');
[prbIdx, dtMt] = obj.pPrbDet(fineCompSignal); 
%% Section 1
for sorter = 1:(length(prbIdx)-1)
 decoder_bin(sorter,:) = fineCompSignal(prbIdx(sorter)+1:prbIdx(sorter+1));
 decoder_message_bin(sorter,:) = fineCompSignal((prbIdx(sorter)+1):(prbIdx(sorter+1)-length(obj.pModulatedHeader)));
 pilot_bin(sorter,:) = fineCompSignal((prbIdx(sorter+1)+1-length(obj.pModulatedHeader)):(prbIdx(sorter+1)));
end
pilot_bin = pilot_bin';
decoder_message_bin=decoder_message_bin';
decoder_bin=decoder_bin';
%%
phaseEst = round(angle(mean(conj(obj.pModulatedHeader).*pilot_bin(:,1)))*2/pi)/2*pi;
phShiftedData = decoder_message_bin(:,1).* exp(-1i*phaseEst);
obj.pQPSKDemodulator = comm.QPSKDemodulator('PhaseOffset', pi/4, ...
 'BitOutput', true);
demodOut_2 = obj.pQPSKDemodulator(phShiftedData);
obj.pDescrambler = comm.Descrambler(prmQPSKReceiver.ScramblerBase, ...
 prmQPSKReceiver.ScramblerPolynomial, prmQPSKReceiver.ScramblerInitialConditions);
deScrData_2 = obj.pDescrambler(demodOut_2);


BER = (length(message_bits)-sum(deScrData_2==message_bits))/(length(message_bits));

received_char = dec2bin(deScrData_2);
received_char = received_char';
symbol8bits = [];
counter = 0;
for bit = 1:length(received_char)/8
 symbol8bits = [symbol8bits; received_char(bit+counter:bit+counter+7)];
 counter = counter + 7;
end
%%
decimal_decoded_data = bin2dec(symbol8bits);
%%
size(decimal_decoded_data);
size(index);
quantizedError = decimal_decoded_data - index;
%%
find(quantizedError~=0);
%%
receivedQuants = codeBook(decimal_decoded_data + 1);
%%
maxofreceivedQuants = max(receivedQuants);
expandedSignal = compand(receivedQuants,Mu,maxofreceivedQuants,'mu/expander');
%%
%Recover the original signal by resampling the expanded signal to the correct value
%Play the recording
samplingFrequency = 16000;
recoveredSignalToTest = resample(expandedSignal,3,2);
sound(recoveredSignalToTest,samplingFrequency)
pause(5)