
load('myfile.mat');
load('OFDM.mat');
%%
length_header = 20;
baseband_sampling = 60e6;
symbol_rate = 2e6;
length_symbol = length(stream);

SimParams.FrameTime = length(Transmitted_signal)./baseband_sampling;
SimParams.StopTime  = 100;


%%
SimParams.Address = 'usb:0'

radio = sdrtx('Pluto');
    radio.RadioID               = 'usb:0';
    radio.CenterFrequency       = 2.2e9;
    radio.BasebandSampleRate    = 60e6;
    radio.SamplesPerFrame       = length(Transmitted_signal);
    radio.Gain                  = 0;
    %%
currentTime = 0;
disp('Transmission has started') % nothing to change here
    % Transmission Process
while currentTime < 100
    data = conj(Transmitted_signal');
    % Data transmission
    step(radio, data);
%disp('sent data')
    % Update simulation time 
    currentTime = currentTime + SimParams.FrameTime;
end

if currentTime ~= 0
    disp('Transmission has ended') % nothing to change here
end    
release(radio);

