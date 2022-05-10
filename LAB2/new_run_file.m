function recorded_signal = new_run_file(prmQPSKReceiver,printData)

%   Copyright 2017 The MathWorks, Inc.
persistent rx radio;
if isempty(rx)
    rx  = QPSKReceiver(...
        'ModulationOrder',                      prmQPSKReceiver.ModulationOrder, ...
        'SampleRate',                           prmQPSKReceiver.Fs, ...
        'DecimationFactor',                     prmQPSKReceiver.Decimation, ...
        'FrameSize',                            prmQPSKReceiver.FrameSize, ...
        'HeaderLength',                         prmQPSKReceiver.HeaderLength, ...
        'NumberOfMessage',                      prmQPSKReceiver.NumberOfMessage, ...
        'PayloadLength',                        prmQPSKReceiver.PayloadLength, ...
        'DesiredPower',                         prmQPSKReceiver.DesiredPower, ...
        'AveragingLength',                      prmQPSKReceiver.AveragingLength, ...
        'MaxPowerGain',                         prmQPSKReceiver.MaxPowerGain, ...
        'RolloffFactor',                        prmQPSKReceiver.RolloffFactor, ...
        'RaisedCosineFilterSpan',               prmQPSKReceiver.RaisedCosineFilterSpan, ...
        'InputSamplesPerSymbol',                prmQPSKReceiver.Interpolation, ...
        'MaximumFrequencyOffset',               prmQPSKReceiver.MaximumFrequencyOffset, ...
        'PostFilterOversampling',               prmQPSKReceiver.Interpolation/prmQPSKReceiver.Decimation, ...
        'PhaseRecoveryLoopBandwidth',           prmQPSKReceiver.PhaseRecoveryLoopBandwidth, ...
        'PhaseRecoveryDampingFactor',           prmQPSKReceiver.PhaseRecoveryDampingFactor, ...
        'TimingRecoveryDampingFactor',          prmQPSKReceiver.TimingRecoveryDampingFactor, ...
        'TimingRecoveryLoopBandwidth',          prmQPSKReceiver.TimingRecoveryLoopBandwidth, ...
        'TimingErrorDetectorGain',              prmQPSKReceiver.TimingErrorDetectorGain, ...
        'PreambleDetectorThreshold',            prmQPSKReceiver.PreambleDetectorThreshold, ...
        'DescramblerBase',                      prmQPSKReceiver.ScramblerBase, ...
        'DescramblerPolynomial',                prmQPSKReceiver.ScramblerPolynomial, ...
        'DescramblerInitialConditions',         prmQPSKReceiver.ScramblerInitialConditions,...
        'BerMask',                              prmQPSKReceiver.BerMask,...)
        'PrintOption',                          printData);
    
    % Create and configure the Pluto System object.
    radio = sdrrx('Pluto');
    radio.RadioID               = 'usb:0';
    radio.CenterFrequency       = prmQPSKReceiver.PlutoCenterFrequency;
    radio.BasebandSampleRate    = prmQPSKReceiver.PlutoFrontEndSampleRate;
    radio.SamplesPerFrame       = prmQPSKReceiver.PlutoFrameLength;
    radio.GainSource            = 'Manual';
    radio.Gain                  = prmQPSKReceiver.PlutoGain;
    radio.OutputDataType        = 'double';
end

% Initialize variables
currentTime = 0;
recorded_signal = [];
rcvdSignal = complex(zeros(prmQPSKReceiver.PlutoFrameLength,1));
while currentTime < prmQPSKReceiver.StopTime
 % Receive signal from the radio
 rcvdSignal = radio();

 % Decode the received message
 %[~,~,~, BER] = rx(rcvdSignal);
recorded_signal = [recorded_signal,rcvdSignal];
 % Update simulation time
 currentTime=currentTime+(radio.SamplesPerFrame / radio.BasebandSampleRate);
end

release(rx);
release(radio);