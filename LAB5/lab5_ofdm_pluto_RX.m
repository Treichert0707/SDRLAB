
radio = sdrrx('Pluto')
    radio.RadioID               = 'usb:0'
    radio.CenterFrequency       = 2.2e9-3e4;
    radio.BasebandSampleRate    = 60e6
    radio.SamplesPerFrame       = length(Transmitted_signal);
    radio.GainSource            = 'Manual'
    radio.Gain                  = 60; % for now...
    radio.OutputDataType        = 'double' % don't change this! 

% Initialize variables
currentTime = 0;
rcvdSignal = complex(zeros(radio.SamplesPerFrame,1));
dummy = 1;
while currentTime < 0.01
    % Receive signal from the radio
    rcvdSignal = radio();
    % Decode the received message
   Received_before_AGC(:,dummy) = rcvdSignal;
   dummy = dummy + 1
    % Update simulation time
    currentTime=currentTime+(radio.SamplesPerFrame / radio.BasebandSampleRate)
end
release(radio);

%%
Received_before_AGC = reshape(Received_before_AGC,[],1);
Received_before_AGC = Received_before_AGC(1:10*length(Transmitted_signal)); % chosing to store 10 samples at max

obj.pAGC = comm.AGC( ...
                'DesiredOutputPower',       2, ...
                'AveragingLength',          50, ...
                'MaxPowerGain',             60);

Received_before_AGC = Received_before_AGC(:);
bufferSignal = Received_before_AGC;
OFDM_received_signal = obj.pAGC(bufferSignal);
%to be used in the receiver...
%% 
